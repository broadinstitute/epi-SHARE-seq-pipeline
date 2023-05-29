library(reticulate)

read_h5ad <- function(
  filename,
  backed = NULL
) {
  python_anndata <- reticulate::import("anndata", convert = FALSE)
  filename <- normalizePath(filename, mustWork = FALSE)
  py_to_r_ifneedbe(python_anndata$read_h5ad(
    filename = filename,
    backed = backed
  ))
}

py_to_r_ifneedbe <- function(x) {
  if (inherits(x, "python.builtin.object")) {
    py_to_r(x)
  } else {
    x
  }
}

#' @name r-py-conversion
#' @export
py_to_r.pandas.core.indexes.base.Index <- function(x) {
  python_builtins <- reticulate::import_builtins()
  out <- python_builtins$list(x)
  attr(out, "name") <- py_to_r_ifneedbe(x$name)
  out
}

#' Convert between Python and R objects
#'
#' @param x A Python object.
#' @param name A name
#' @param value A value
#'
#' @return An \R object, as converted from the Python object.
#'
#' @name r-py-conversion
#' @export
`[[<-.collections.abc.MutableMapping` <- function(x, name, value) {
  if (!is.null(value)) {
    reticulate::py_set_item(x, name, value)
  } else if (name %in% x$keys()) {
    reticulate::py_del_item(x, name)
  }
}

#' @name r-py-conversion
#' @export
`[[.collections.abc.Mapping` <- function(x, name) {
  if (name %in% x$keys()) {
    py_to_r_ifneedbe(reticulate::py_get_item(x, name))
  } else {
    NULL
  }
}

#' @name r-py-conversion
#' @export
`[<-.collections.abc.MutableMapping` <- `[[<-.collections.abc.MutableMapping`
#
#' @name r-py-conversion
#' @export
`[.collections.abc.Mapping` <- `[[.collections.abc.Mapping`
#
#' @name r-py-conversion
#' @export
`names.collections.abc.Mapping` <- function(x) {
  python_builtins <- reticulate::import_builtins()
  python_builtins$list(x$keys())
}

#' @name r-py-conversion
#' @export
`py_to_r.collections.abc.Set` <- function(x) {
  python_builtins <- reticulate::import_builtins()
  python_builtins$list(x)
}

#' @name r-py-conversion
#' @export
py_to_r.pandas.core.indexes.base.Index <- function(x) {
  python_builtins <- reticulate::import_builtins()
  out <- python_builtins$list(x)
  attr(out, "name") <- py_to_r_ifneedbe(x$name)
  out
}

#' @name r-py-conversion
#' @export
py_to_r.collections.abc.KeysView <- function(x) {
  python_builtins <- reticulate::import_builtins()
  python_builtins$list(x)
}

#' @name r-py-conversion
#' @export
`py_to_r.collections.abc.Mapping` <- function(x) {
  python_builtins <- reticulate::import_builtins()

  x_list <- python_builtins$dict(x)

  # convert members of x_list if need be
  for (i in seq_along(x_list)) {
    if (inherits(x_list[[i]], "python.builtin.object")) {
      x_list[[i]] <- py_to_r_ifneedbe(x_list[[i]])
    }
  }

  x_list
}


#' @importFrom Matrix sparseMatrix
py_to_r.scipy.sparse.csc.csc_matrix <- function(x) {
  Matrix::sparseMatrix(
    i = as.integer(py_to_r_ifneedbe(x$indices))+1,
    p = as.integer(py_to_r_ifneedbe(x$indptr)),
    x = as.vector(py_to_r_ifneedbe(x$data)),
    dims = as.integer(dim(x))
  )
}