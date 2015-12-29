#' Thresholded distances between columns of a matrix.
#'
#' Unlike the \code{\link{dist}} function which computes
#' distances between matrix _rows_, \code{tdist} computes distances
#' exceeding a threshold between matrix _columns_.
#'
#' Increase p to cut down the total number of candidate pairs evaluated,
#' at the expense of costlier truncated SVDs.
#'
#' This function doesn't work as well generally yet as the \code{\link{tcor}}
#' companion function in this package.
#'
#' @param A an m by n real-valued dense or sparse matrix
#' @param t a threshold distance value, if missing an estimate derived from a 1-d SVD projection will be used
#' @param p projected subspace dimension
#' @param filter "local" filters candidate set sequentially,
#'  "distributed" computes thresholded correlations in a parallel code section which can be
#'  faster but requires that the data matrix is available (see notes).
#' @param method the distance measure to be used, one of
#'          "euclidean", or "manhattan".
#' Any unambiguous substring can be given.
#' @param dry_run set \code{TRUE} to return statistics and truncated SVD for tuning
#' \code{p} (see notes)
#' @param restart either output from a previous run of \code{tdist} with \code{dry_run=TRUE},
#' or direct output from from \code{\link{irlba}} used to restart the \code{irlba}
#' algorithm when tuning \code{p} (see notes)
#' @param ... additional arguments passed to \code{\link{irlba}}
#'
#' @return A list with elements:
#' \enumerate{
#'   \item \code{indices} A three-column matrix. The  first two columns contain
#'         indices of vectors meeting the distance threshold \code{t},
#'         the third column contains the corresponding distance value.
#'   \item \code{longest_run} The largest number of successive entries in the
#'     ordered first singular vector within a projected distance defined by the
#'     correlation threshold.
#'   \item \code{n} The total number of _possible_ vectors that meet
#'     the correlation threshold identified by the algorithm.
#'   \item \code{total_time} Total run time.
#' }
#' @importFrom irlba irlba
#' @importFrom stats dist
#' @export
tdist = function(A, t, p=10,
                 filter=c("distributed", "local"),
                 method=c("euclidean", "manhattan"),
                 dry_run=FALSE, restart, t_est, ...)
{
  filter = match.arg(filter)
  method = match.arg(method)
  if(ncol(A) < p) p = max(1, floor(ncol(A) / 2 - 1))
  t0 = proc.time()
  if(missing(restart)) L  = irlba(A, p, ...)
  else
  {
    # Handle either output from tcor(..., dry_run=TRUE), or direct output from irlba:
    if("restart" %in% names(restart)) restart = restart$restart
    L = irlba(A, p, v=restart, ...)
  }
  t1 = (proc.time() - t0)[[3]]
  if(missing(t))
  {
    # Estimate a threshold based on a 1-d projection
# XXX improve me
    v = L$v[order(L$v[,1]), 1]
    if(missing(t_est)) t_est = 0.01
    t = quantile(L$d[1] * (v - v[1]), probs=t_est)
  }

  normlim = switch(method,
                   euclidean = t ^ 2,
                   maximum = nrow(A) * t ^ 2,  # XXX unlikely to be a good bound improve? XXX
                   manhattan   = t ^ 2)        # just bound by 2-norm
  full_dist_fun =
    switch(method,
           euclidean = function(idx) vapply(1:nrow(idx), function(k) sqrt(crossprod(A[, idx[k,1]] - A[, idx[k, 2]])), 1),
           manhattan = function(idx) vapply(1:nrow(idx), function(k) sum(abs(A[, idx[k,1]] - A[, idx[k, 2]])), 1),
           maximum = function(idx) vapply(1:nrow(idx), function(k) max(abs(A[, idx[k,1]] - A[, idx[k, 2]])), 1)
  )
  filter_fun =  function(v, t) v <= t

  ans = two_seven(A, L, t, filter, normlim=normlim, full_dist_fun=full_dist_fun, filter_fun=filter_fun, dry_run=dry_run)
  if(dry_run) return(list(restart=L, longest_run=ans$longest_run, n=ans$n, t=t, svd_time=t1))
  return(c(ans, svd_time=t1, total_time=(proc.time() - t0)[[3]]))
}
