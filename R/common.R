#' linear time longest run search (A. Poliakov), find the longest
#' run of values in the vector within the specified distance
#' @param v a vector with entries ordered in increasing order
#' @param limit distance interval
#' @return run length
#' @keywords internal
longrun = function(v, limit)
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
      while(lower < upper && v[upper] - v[lower] > limit) lower = lower + 1
    }
  }
  ell
}

#' Steps 2--7 of Algorithm 2.1, factored into a common function that can be used by a variety of distance metrics
#' @param A data matrix
#' @param L truncated SVD of A
#' @param t scalar threshold value
#' @param filter "distributed" for full threshold evaluation of pruned set on parallel workers, "local" for sequential evaluation of full threshold of pruned set to avoid copying data matrix.
#' @param normlim the squared norm limit in step 4, default value is for correlation
#' @param full_dist_fun non-projected distance function of a two-column matrix of rows of column indices that needs scoped access to A (step 7), default function is for correlation
#' @param filter_fun filter function of a vector and scalar that thresholds vector values from full_dist_fun, returning a logical vector of same length as v (step 7), default function is for correlation
#' @param dry_run a logical value, if \code{TRUE} quickly return statistics useful for tuning \code{p}
#' @return a list with indices, ell, n, and longest_run entries, unless dry_run=\code{TRUE} in which case
#' a list with ell and n is returned
#' @importFrom foreach foreach %dopar%
#' @keywords internal
two_seven = function(A, L, t, filter=c("distributed", "local"), normlim=2 * (1 - t),
                     full_dist_fun=function(idx) vapply(1:nrow(idx), function(k) cor(A[, idx[k,1]], A[, idx[k, 2]]), 1),
                     filter_fun=function(v, t) v >= t, dry_run=FALSE)
{
  filter = match.arg(filter)
# Find the projection among the first few basis vectors with the shortest
# maximum run length to minimize work in the next step. This is a cheap but
# usually not very significant optimization.
  p = length(L$d)
  ells = lapply(1:min(2,p), function(N)
  {
    P = order(L$v[, N])
    limit = sqrt(normlim) / L$d[N]
    ell = longrun(L$v[order(L$v[, N]), N], limit)
    list(P=P, limit=limit, ell=ell, N=N)
  })
  ellmin = which.min(vapply(ells, function(x) x$ell, 1))
  P = ells[[ellmin]]$P
  limit = ells[[ellmin]]$limit
  ell = min(ells[[ellmin]]$ell, ncol(A) - 1)

  if(dry_run)
  {
    d = diff(L$v[P, 1:p, drop=FALSE], lag=1) ^ 2 %*% L$d[1:p] ^ 2
    return(list(longest_run=ell, n=sum(d <= normlim)))
  }

# The big union in step 4 of algorithm 2.1 follows, combined with step 6 to
# convert back to original indices, and step 7 to evaluate the candiadtes.
# Each step from 1 to ell is independent of the others; the steps can run
# in parallel.
  combine = function(x, y)
  {
    list(indices=rbind(x$indices, y$indices), n=x$n + y$n)
  }

  if(filter == "distributed")
  {
    indices = foreach(i=1:ell, .combine=combine, .inorder=FALSE) %dopar%
    {
      d = diff(L$v[P, 1:p, drop=FALSE], lag=i) ^ 2 %*% L$d[1:p] ^ 2
      # These ordered indices meet the projected threshold:
      j = which(d <= normlim)
      n = length(j)
      # return original un-permuted column indices that meet true threshold
      # (step 7), including the number of possible candidates for info.:
      if(n == 0)
      {
        ans = vector("list", 2)
        names(ans) = c("idx", "n")
        ans$n = n
        return(ans)
      }
      v = full_dist_fun(cbind(P[j], P[j + i]))
      h = filter_fun(v, t)
      j = j[h]
      v = v[h]
      return(list(indices=cbind(i=P[j], j=P[j + i], val=v), n=n))
    }
    return(c(indices, longest_run=ell))
  }

# filter == "local" case, preventing copy of the data matrix to the workers
  indices = foreach(i=1:ell, .combine=combine, .inorder=FALSE, .noexport="A") %dopar%
  {
    d = diff(L$v[P, 1:p, drop=FALSE], lag=i) ^ 2 %*% L$d[1:p] ^ 2
    # These ordered indices meet the projected threshold:
    j = which(d <= normlim)
    n = length(j)
    # return original un-permuted column indices that meet true threshold
    # (step 7), including the number of possible candidates for info.:
    if(n == 0)
    {
      ans = vector("list", 2)
      names(ans) = c("indices", "n")
      ans$n = n
      return(ans)
    }
    list(indices=cbind(i=P[j], j=P[j + i]), n=n)
  }
  v = full_dist_fun(indices$indices)
  h = filter_fun(v, t)
  indices$indices = cbind(indices$indices[h,], val=v[h])
  c(indices, longest_run=ell)
}
