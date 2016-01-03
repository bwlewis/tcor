#' Compute thresholded distances between rows or columns of a matrix
#'
#' Compute and return distances and indices of rows or columns within a specified distance threshold
#' with respect to a specified distance metric. The algorithm works best for Euclidean
#' distance (the default option).
#' Alternatively compute the \code{t} closest rows when \code{rank=TRUE}. Or use
#' \code{columns=TRUE} to compute distances between columns instead, which is somewhat
#' cheaper for this algorithm than computing row distances.
#' Increase p to cut down the total number of candidate pairs evaluated,
#' at the expense of costlier truncated SVDs.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold distance value either in absolute distance (the default) or rank order (see \code{rank} below);
#' if missing an estimate derived from a 1-d SVD projection will be used
#' @param p projected subspace dimension
#' @param filter "local" filters candidate set sequentially,
#'  "distributed" computes thresholded correlations in a parallel code section which can be
#'  faster but requires that the data matrix is available (see notes).
#' @param method the distance measure to be used, one of
#'          "euclidean", or "manhattan".
#' Any unambiguous substring can be given.
#' @param rank when \code{TRUE}, the threshold \code{t} represents the top \code{t}
#' closest vectors, otherwise the threshold \code{t} specifies absolute distance; when
#' \code{rank=TRUE} then \code{t} must also be specified
#' @param dry_run set \code{TRUE} to return statistics and truncated SVD for tuning
#' \code{p} (see notes)
#' @param max_iter when \code{rank=TRUE}, a portion of the algorithm may iterate; this
#' number sets the maximum numer of such iterations
#' @param columns set to \code{TRUE} to compute distances between matrix columns instead
#' of rows, saving the expense of a matrix transpose (which can be significant if \code{A} is large)
#' @param restart either output from a previous run of \code{tdist} with \code{dry_run=TRUE},
#' or direct output from from \code{\link{irlba}} used to restart the \code{irlba}
#' algorithm when tuning \code{p} (see notes)
#' @param ... additional arguments passed to \code{\link{irlba}}
#'
#' @return A list with elements:
#' \enumerate{
#'   \item \code{indices} A three-column matrix. The  first two columns contain
#'         indices of rows meeting the distance threshold \code{t},
#'         the third column contains the corresponding distance value (not returned
#'         when \code{dry_run=TRUE}).
#'   \item \code{restart} A truncated SVD returned by the IRLBA used to restart the
#'   algorithm (only returned when \code{dry_run=TRUE}).
#'   \item \code{n} The total number of _possible_ vectors that meet
#'     the correlation threshold identified by the algorithm.
#'   \item \code{longest_run} The largest number of successive entries in the
#'     ordered first singular vector within a projected distance defined by the
#'     correlation threshold; Equivalently, the number of \code{n * p} matrix
#'     vector products employed in the algorithm, not counting the truncated SVD step.
#'   \item \code{t} The threshold value.
#'   \item \code{svd_time} Time to compute truncated SVD.
#'   \item \code{total_time} Total run time.
#' }
#'
#' @note When \code{rank=TRUE} the method returns at least, and perhaps more than, the top \code{t} closest
#' indices and their distances, unless they could not be found within the iteration
#' limit \code{max_iter}.
#' @seealso \code{\link{dist}}, \code{\link{tcor}}
#' @references \url{http://arxiv.org/abs/1512.07246} (preprint)
#' @examples
#' x <- matrix(rnorm(100*20), nrow=100)
#' # Find the top 10 closest vectors with respect to Euclidean distance:
#' td <- tdist(x, 10, rank=TRUE)
#' print(td$indices[1:10,])
#'
#' # Compare with distances from `dist`:
#' d <- dist(x)
#' print(sort(d)[1:10])
#'
#' @importFrom irlba irlba
#' @importFrom stats dist
#' @export
tdist = function(A, t, p=10,
                 filter=c("distributed", "local"),
                 method=c("euclidean", "manhattan", "maximum"), rank=FALSE,
                 dry_run=FALSE, max_iter=4, columns=FALSE, restart, ...)
{
  filter = match.arg(filter)
  method = match.arg(method)
  if(!columns) A = base::t(A)  # XXX expensive, find a better approach...
  if(ncol(A) < p) p = max(1, floor(ncol(A) / 2 - 1))

  nlim = function(t)
           switch(method,
                  euclidean = t ^ 2,
                  maximum = nrow(A) * t ^ 2,  # XXX unlikely to be a good bound?
                  manhattan   = t ^ 2)        # just bound by 2-norm, not so great either

  t0 = proc.time()
  if(missing(restart)) L  = irlba(A, p, ...)
  else
  {
    # Handle either output from tcor(..., dry_run=TRUE), or direct output from irlba:
    if("restart" %in% names(restart)) restart = restart$restart
    L = irlba(A, p, v=restart, ...)
  }
  t1 = (proc.time() - t0)[[3]]
  if(missing(t) && rank) stop("t must be specified when rank=TRUE")
  N = 1
  if(rank) N = t
  if(missing(t) || rank)
  {
    # Estimate a threshold based on a 1-d projection, with crude tuning over a range of values.
    v = L$v[order(L$v[,1]), 1]
    ts = sort(c(quantile(L$d[1] * (v[-1] - v[1]), probs=c(0.001, 0.01, 0.1)),
              quantile(L$d[1] * (v[length(v)] - v[-length(v)]), probs=c(0.001, 0.01, 0.1))))
    # (that last expression considers two possible SVD bases of different sign)
    as = lapply(ts, function(t) two_seven(A,L,t,filter,normlim=nlim(t),dry_run=TRUE)$n)
    i = which(as > 0)
    if(length(i) > 0) t = ts[min(i)]
    else t = ts[3]
    attr(t, "names") = c()
    if(rank && method == "maximum") t = t / 4    # just a fudge factor, these bounds are not great
    if(rank && method == "manhattan") t = t * 4
  }

  full_dist_fun =
    switch(method,
           euclidean = function(idx) vapply(1:nrow(idx), function(k) sqrt(crossprod(A[, idx[k,1]] - A[, idx[k, 2]])), 1),
           manhattan = function(idx) vapply(1:nrow(idx), function(k) sum(abs(A[, idx[k,1]] - A[, idx[k, 2]])), 1),
           maximum = function(idx) vapply(1:nrow(idx), function(k) max(abs(A[, idx[k,1]] - A[, idx[k, 2]])), 1)
  )
  filter_fun =  function(v, t) v <= t

  iter = 1
  while(iter <= max_iter)
  {
    ans = two_seven(A, L, t, filter, normlim=nlim(t), full_dist_fun=full_dist_fun, filter_fun=filter_fun, dry_run=dry_run)
    if(dry_run) return(list(restart=L, longest_run=ans$longest_run, n=ans$n, t=t, svd_time=t1))
    if(!rank || (nrow(ans$indices) >= N)) break
    iter = iter + 1
# back off faster as we get closer to avoid too much filtering, at the expense of maybe more iterations
    t = t * (2 - (nrow(ans$indices)/N)^(1/4)) 
  }
  ans$indices = ans$indices[order(ans$indices[,"val"]),]
  c(ans, svd_time=t1, total_time=(proc.time() - t0)[[3]])
}
