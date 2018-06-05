#' @include derFunc.R
# @include derringerLs.R
# @include derringerRs.R
# @include derringerTs.R

NULL
###############################################
# --------------------------------------------#
# Class desFunc                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Class definition for desFunc
# --------------------------------------------

# desFunc class
# 
# @name desFunc
setClassUnion("desFunc", c("derFunc"))

# --------------------------------------------
# Accesssor functions for desFunc
# --------------------------------------------

#' Method defining the $ operator for the desFunc class
#' 
#' @inheritParams overview
setMethod("$", "desFunc",
          function(x, name) slot(x, name))

