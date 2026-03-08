###############################################################################
# zzz.R
# Package startup hooks
# ---------------------------------------------------------------------------
# Registers S3 methods for generics defined in Suggests packages
# (ggplot2, generics) so that dispatch works when those packages are loaded.
###############################################################################

.onLoad <- function(libname, pkgname) {
  # Register autoplot.svyder with ggplot2::autoplot
  .s3_register("ggplot2::autoplot", "svyder")

  # Register tidy/glance with generics package
  .s3_register("generics::tidy", "svyder")
  .s3_register("generics::glance", "svyder")

  invisible()
}


# Simplified s3_register, adapted from vctrs::s3_register().
# Conditionally registers an S3 method for a generic in a Suggests package.
.s3_register <- function(generic, class, method = NULL) {
  pieces <- strsplit(generic, "::")[[1]]
  package <- pieces[[1]]
  generic_name <- pieces[[2]]

  caller_env <- parent.frame()

  get_method <- function() {
    if (!is.null(method)) return(method)
    method_name <- paste0(generic_name, ".", class)
    top_env <- topenv(caller_env)
    if (exists(method_name, envir = top_env, inherits = FALSE)) {
      get(method_name, envir = top_env)
    } else {
      NULL
    }
  }

  register <- function(...) {
    if (!isNamespaceLoaded(package)) return(invisible())
    envir <- asNamespace(package)
    method_fn <- get_method()
    if (is.null(method_fn)) return(invisible())
    if (exists(generic_name, envir = envir)) {
      registerS3method(generic_name, class, method_fn, envir = envir)
    }
    invisible()
  }

  # Register immediately if package already loaded, otherwise defer
  if (isNamespaceLoaded(package)) {
    register()
  } else {
    setHook(
      packageEvent(package, "onLoad"),
      function(...) register()
    )
  }

  invisible()
}
