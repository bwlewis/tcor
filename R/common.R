#' linear time longest run search: find the longest run of values in the
#' ordered vector within the specified limit
#' @param v a vector with entries ordered in increasing order
#' @param limit distance interval
#' @param group optional vector with entries -1 and 1 corresponding to the group membership of each element in \code{v} (two groups)
#' @return run length
#' @keywords internal
longrun = function(v, limit, group=NULL)
{
  lower <- 1
  ell = 1
  for(upper in 2:length(v))
  {
    if(v[upper] - v[lower] <= limit)
    {
      ell = max(ell, upper - lower + 1)
    } else
    {
      if(is.null(group)) while (lower < upper && v[upper] - v[lower] > limit) lower = lower + 1
      else while (lower < upper && v[upper] - v[lower] > limit && prod(group[lower:upper]) < 0) lower = lower + 1
    }
  }
  ell
}

#' Steps 2--7 of Algorithm 2.1, factored into a common function that can be used by a variety of distance metrics
#' @param A data matrix
#' @param L truncated SVD of A
#' @param t scalar threshold value
#' @param filter "distributed" for full threshold evaluation of pruned set on parallel workers,
#'         "local" for sequential evaluation of full threshold of pruned set to avoid copying data matrix.
#' @param normlim the squared norm limit in step 4, default value is for correlation
#' @param full_dist_fun non-projected distance function of a two-column matrix of rows of column
#'        indices that needs scoped access to A (step 7), default function is for correlation
#' @param filter_fun filter function of a vector and scalar that thresholds vector values
#'        from full_dist_fun, returning a logical vector of same length as v (step 7), default function is for correlation
#' @param dry_run a logical value, if \code{TRUE} quickly return statistics useful for tuning \code{p}
#' @param anti a logical value, if \code{TRUE} also include anti-correlated vectors
#' @param group either \code{NULL} for no grouping, or a vector of length \code{ncol(A)} consisting of \code{-1, 1} values
#'        indicating group membership of the columns.
#' @return a list with indices, ell, tot, and longest_run entries, unless dry_run=\code{TRUE} in which case
#' a list with ell and tot is returned
#' @importFrom foreach foreach %dopar% %do%
#' @keywords internal
two_seven = function(A, L, t, filter=c("distributed", "local"), normlim=2 * (1 - t),
                     full_dist_fun=function(idx) vapply(1:nrow(idx), function(k) cor(A[, idx[k,1]], A[, idx[k, 2]]), 1),
                     filter_fun=function(v, t) v >= t, dry_run=FALSE, anti=FALSE, group=NULL)
{
  filter = match.arg(filter)
  nx = ncol(A)
  if(!is.null(group)) nx = sum(group == group[1])  # number of columns in 1st array
  grouped = ! is.null(group)
# Find the projection among the first few basis vectors with the shortest
# maximum run length to minimize work in the next step. This is a cheap but
# usually not very significant optimization.
  p = length(L$d)
  ells = lapply(1:min(2, p), function(N)
  {
    P = order(L$v[, N])
    limit = sqrt(normlim) / L$d[N]
    ell = longrun(L$v[order(L$v[, N]), N], limit, group[P])
    list(P=P, limit=limit, ell=ell, N=N)
  })
  ellmin = which.min(vapply(ells, function(x) x$ell, 1))
  P = ells[[ellmin]]$P
  ell = min(ells[[ellmin]]$ell, ncol(A) - 1)
  eN = ells[[ellmin]]$N
  elim = ells[[ellmin]]$limit

# In the include anticorrelated case, we can use the same permutation but
# likely increase ell.
  if(anti)
  {
    v2 = c(L$v[, eN], -L$v[, eN])
    v2p = order(v2)
    v2 = v2[v2p]
    ell = longrun(v2, sqrt(normlim) / L$d[eN], c(group, group)[v2p])
  }

  if(dry_run)
  {
    d = diff(L$v[P, 1:p, drop=FALSE], lag=1) ^ 2 %*% L$d[1:p] ^ 2
    return(list(longest_run=ell, tot=sum(d <= normlim), t=t))
  }

# The big union in step 4 of algorithm 2.1 follows, combined with step 6 to
# convert back to original indices, and step 7 to evaluate the candiadtes.
# Each step from 1 to ell is independent of the others; the steps can run
# in parallel.
  combine = function(x, y)
  {
    list(indices=rbind(x$indices, y$indices), tot=x$n + y$n)
  }

# codetools has trouble detecting the foreach variable i below. We define
# it here to supress CRAN NOTEs and lintr warnings (is this really a
# problem in foreach or codetools?).
  i = 1

  if(filter == "distributed")
  {
    indices = foreach(i=1:ell, .combine=combine, .inorder=FALSE, .packages=c("tcor", "Matrix")) %dopar%
    {
      d2 = Inf
      # restrict focus to candidates from each group (if specified)
      if(grouped)
      {
        # (all this does is make the matrix vector product cheaper)
        gidx = which(diff(group[P], lag=i) != 0)
        D = L$v[P[gidx + i], 1:p] - L$v[P[gidx], 1:p]
        if(anti) D2 = L$v[P[gidx + i], 1:p] + L$v[P[gidx], 1:p]
      } else
      {
        D = diff(L$v[P, 1:p, drop=FALSE], lag=i)
        if(anti) D2 = diffint(L$v[P, 1:p, drop=FALSE], lag=i, sign=1)
      }
      d = D ^ 2 %*% L$d[1:p] ^ 2
      if(anti) d2 = D2 ^ 2 %*% L$d[1:p] ^ 2
      # These ordered indices meet the projected threshold:
      j = c(which(d <= normlim), which(d2 <= normlim))
      if(grouped) j = gidx[j]
      n = length(j)
      # return original un-permuted column indices that meet true threshold
      # (step 7), and the number of possible candidates:
      if(n == 0)
      {
        ans = vector("list", 2)
        names(ans) = c("indices", "tot")
        ans$indices = cbind(i=integer(0), j=integer(0), val=double(0))
        ans$tot = n
        return(ans)
      }
      v = full_dist_fun(cbind(P[j], P[j + i]))
      h = filter_fun(v, t)
      j = j[h]
      v = v[h]
      ans_i = P[j]
      ans_j = P[j + i]
      if(grouped)
      {
        idx = which(ans_i  > nx)
        if(length(idx) > 0)
        {
          j2 = ans_i[idx]
          ans_i[idx] = ans_j[idx]
          ans_j[idx] = j2
        }
        ans_i = ans_i %% nx
        ans_j = ans_j %% nx
        ans_i[ans_i == 0] = nx
        ans_j[ans_j == 0] = nx
      }
      list(indices=cbind(i=ans_i , j=ans_j, val=v), tot=n)
    }
    return(c(indices, longest_run=ell, t=t))
  }

# The filter == "local" case, preventing copy of the data matrix to the workers
  indices = foreach(i=1:ell, .combine=combine, .inorder=FALSE, .noexport="A") %do% # XXX
  {
    d2 = Inf
    # restrict focus to candidates from each group (if specified)
    if(grouped)
    {
      # (all this does is make the matrix vector product cheaper)
      gidx = which(diff(group[P], lag=i) != 0)
      D = L$v[P[gidx + i], 1:p] - L$v[P[gidx], 1:p]
      if(anti) D2 = L$v[P[gidx + i], 1:p] + L$v[P[gidx], 1:p]
    } else
    {
      D = diff(L$v[P, 1:p, drop=FALSE], lag=i)
      if(anti) D2 = diffint(L$v[P, 1:p, drop=FALSE], lag=i, sign=1)
    }
    d = D ^ 2 %*% L$d[1:p] ^ 2
    if(anti) d2 = D2 ^ 2 %*% L$d[1:p] ^ 2
    # These ordered indices meet the projected threshold:
    j = c(which(d <= normlim), which(d2 <= normlim))
    if(grouped) j = gidx[j]
    n = length(j)
    # return original un-permuted column indices that meet true threshold
    # (step 7), including the number of possible candidates for info.:
    if(n == 0)
    {
      ans = vector("list", 2)
      names(ans) = c("indices", "tot")
      ans$tot = n
      return(ans)
    }
    ans_i = P[j]
    ans_j = P[j + i]
    if(grouped)
    {
      idx = which(ans_i  > nx)
      if(length(idx) > 0)
      {
        j2 = ans_i[idx]
        ans_i[idx] = ans_j[idx]
        ans_j[idx] = j2
      }
    }
    list(indices=cbind(i=ans_i, j=ans_j), tot=n)
  }
  v = full_dist_fun(indices$indices)
  h = filter_fun(v, t)
  indices$indices = cbind(indices$indices[h,], val=v[h])
  if(grouped)
  {
    indices$indices[, 1:2] = indices$indices[, 1:2] %% nx
    indices$indices[indices$indices[, 1] == 0, 1] = nx
    indices$indices[indices$indices[, 2] == 0, 2] = nx
  }
  c(indices, longest_run=ell, t=t)
}

# replacement for diff that also supports sums
diffint = function(x, lag=1, sign=-1)
{
  tail(x, nrow(x) - lag) + sign * head(x, nrow(x) - lag)
}
