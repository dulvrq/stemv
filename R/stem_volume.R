
#' A list of functions for calculating stem volume
#'
#' A list of functions for calculating stem volume of trees based on DBH, H, and species.
#' Main functions are [stemVolume] and [volumeName].
#'
#' @name stem-volume-functions
NULL


#' Get coefficients for calculating stem volume
#'
#' @param dt_stem_i  Dataframe from the predefined table.
#' @param D  Numeric. DBH in cm.
#' @return
#'   A list of coefficients
#' @export
#'
getCoefs <- function(dt_stem_i, D){
  wch_c <- which(dt_stem_i$D_lower <= D & dt_stem_i$D_upper > D)
  list_coefs <- list(
    a = dt_stem_i$coef_a[wch_c],
    b = dt_stem_i$coef_b[wch_c],
    c = dt_stem_i$coef_c[wch_c],
    d = dt_stem_i$coef_d[wch_c]
  )
  return(list_coefs)
}


#' Round, not to even digit
#'
#' Implement usual rounding, not to even digits.
#'
#' @param x Numeric.
#' @param digits integer. digits to round.
#' @return Rounded numeric.
#'
round2 <- function (x, digits = 0) return(floor(x * 10^digits + 0.5) / 10^digits)


#' Calculate stem volume using an equation
#'
#' Calculate stem volume using an equation that corresponds with the equation type.
#' This function can handle numeric vectors with different DBH (Diameter at Breast Height) classes.
#'
#' @param dt_stem_i  Dataframe from the predefined table.
#' @param D  Numeric/vector. DBH in cm.
#' @param H  Numeric/vector. Tree height in m.
#' @param eq_type  Numeric for equation type (currently, 1-4 is defined).
#' @param l_c  A list of parameters. If NULL, parameters will be derived from `dt_stem_i`.
#' @return
#'   Calculated stem volume
#'
#' @importFrom purrr
#'    map_dbl map2 pmap
#' @export
#'
calcStemVolume <- function(dt_stem_i, D, H, eq_type = 1, l_c = NULL){
  ## get parameters that correspond with DBH size ---
  if(is.null(l_c)){
    wch_c <- map_dbl(D, function(x) which(dt_stem_i$D_lower <= x & dt_stem_i$D_upper > x))
    ## if the DBH classes are different, get parameters for each D
    if(all(wch_c == wch_c[1])){
      l_c <- getCoefs(dt_stem_i, D[1])
    } else {
      l_c <- map2(list(dt_stem_i), D, getCoefs)
    }
  }

  ## 1) using single set of parameters ---
  if(length(l_c[[1]]) == 1){
    # eq. 1 ---
    if(eq_type == 1){
      ## calculate stem volume ---
      V <- 10^(l_c[["a"]] + l_c[["b"]]*log10(D) + l_c[["c"]]*log10(H))

      # eq. 2 ---
    } else if (eq_type == 2){
      ## DBH unit conversion if required
      if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      ## calculate stem volume ---
      l_c$b[is.na(l_c$b)] <- 0 # for Kochi akamatsu
      V <- l_c[["a"]]*H + l_c[["b"]]*Dm^2 + l_c[["c"]]*Dm^2*H + l_c[["d"]]

      # eq. 3 ---
    } else if (eq_type == 3){
      ## DBH unit conversion if required
      if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      ## calculate stem volume ---
      V <- l_c[["c"]]*Dm^2*H + l_c[["d"]]

      # eq. 4 ---
    } else if (eq_type == 4){
      ## calculate stem volume ---
      V <- 10^(l_c[["a"]] + l_c[["b"]]*log10(D^2*H))

      # eq. 5 ---
      # special case for kami-yaku tennen-sugi, D around 61cm
    } else if (eq_type == 5){
      ## calculate stem volume ---
      V <- 10^(l_c[["a"]] + l_c[["b"]]*log10(D) + l_c[["c"]]*log10(H)) # same as eq. 1
      if(any(D >= 61)){
        wch_replace <- which(D >= 61)
        V[wch_replace] <- 0.1513 + 0.8898 * V[wch_replace] # coefficients from kami-yaku tennen-sugi
      }

    } else {
      stop("No equation.")
    }

    ## 2) using a set of parameters for each D, a little slower ---
  } else {
    # eq. 1 ---
    if(eq_type == 1){
      V <- pmap(list(D = D, H = H, l_c = l_c),
                function(D, H, l_c) 10^(l_c[["a"]] + l_c[["b"]]*log10(D) + l_c[["c"]]*log10(H))) |>
        unlist()

      # eq. 2 ---
    } else if (eq_type == 2){
      ## DBH unit conversion if required
      if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      ## calculate stem volume ---
      V <- pmap(list(Dm = Dm, H = H, l_c = l_c),
                function(Dm, H, l_c){
                  l_c$b[is.na(l_c$b)] <- 0 # for Kochi akamatsu
                  l_c[["a"]]*H + l_c[["b"]]*Dm^2 + l_c[["c"]]*Dm^2*H + l_c[["d"]]
                }) |>
        unlist()

      # eq. 3 ---
    } else if (eq_type == 3){
      ## DBH unit conversion if required
      if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      ## calculate stem volume ---
      V <- pmap(list(Dm = Dm, H = H, l_c = l_c),
                function(Dm, H, l_c) l_c[["c"]]*Dm^2*H + l_c[["d"]]) |>
        unlist()


      # eq. 4 ---
    } else if (eq_type == 4){
      ## calculate stem volume ---
      V <- pmap(list(D = D, H = H, l_c = l_c),
                function(D, H, l_c) 10^(l_c[["a"]] + l_c[["b"]]*log10(D^2*H))) |>
        unlist()

      # eq. 5 ---
      # special case for kami-yaku tennen-sugi, D around 61cm
    } else if (eq_type == 5){
      ## calculate stem volume ---
      V <- pmap(list(D = D, H = H, l_c = l_c),
                function(D, H, l_c) 10^(l_c[["a"]] + l_c[["b"]]*log10(D) + l_c[["c"]]*log10(H))) |>
        unlist() # same as eq. 1
      if(any(D >= 61)){
        wch_replace <- which(D >= 61)
        V[wch_replace] <- 0.1513 + 0.8898 * V[wch_replace] # coefficients from kami-yaku tennen-sugi
      }

    } else {
      stop("No equation.")
    }
  }

  return(V)
}


#' Linear interpolation
#'
#' Linear imputation between 2-diameter & volume.
#' This can be used not only for DBH & V, but for other parameters.
#'
#' @param d1  Numeric/vector. Lower bound DBH in cm.
#' @param v1  Numeric/vector. Calculated volume that corresponds with `d1`.
#' @param d2  Numeric/vector. Upper bound DBH in cm.
#' @param v2  Numeric/vector. Calculated Volume that corresponds with `d2`.
#' @param d  Numeric/vector. DBH in cm.
#' @return
#'   Linear interpolated stem volume or other variables
#'
linearImpute <- function(d1, v1, d2, v2, d){
  return(v1 + (v2 - v1) / (d2 - d1) * (d - d1))
}

#' Linear interpolation for volume
#'
#' Linear imputation between 2-diameter & volume.
#' Both for Volume calculation & linear interpolation.
#'
#' @param dt_stem_i Dataframe from the predefined table.
#' @param d1  Numeric/vector. Lower bound DBH in cm.
#' @param d2  Numeric/vector. Upper bound DBH in cm.
#' @param D  Numeric/vector. DBH in cm.
#' @param H  Numeric/vector. Tree height in m.
#' @param eq_type  Numeric for equation type (currently, 1-4 is defined).
#' @return
#'   Calculated stem volume
#'
linearImputeVolume <- function(dt_stem_i, d1, d2, D, H, eq_type){
  v1 <- calcStemVolume(dt_stem_i, d1, H, eq_type)
  v2 <- calcStemVolume(dt_stem_i, d2, H, eq_type)
  return(v1 + (v2 - v1) / (d2 - d1) * (D - d1))
}


#' Calculate stem volume based on junction adjustment, for single value.
#'
#' @description
#' Calculate stem volume using an equation corresponding with the equation type.
#' \code{calcStemVolumeAdj()} only accepts single numeric value for D & H.
#' To handle numeric vector, use [calcStemVolumeMulti()].
#'
#' @details
#' This is a wrapper function to apply 3/5- points moving window average for calculating stem volume
#' when DBH is located at the junction of equations.\cr
#' VERY SLOW for the large amount of data.
#'
#' See the following reference for calculation details:\cr
#'  Hosoda et al. (2010):
#' "Differences between the present stem volume tables and the values of the volume equations, and their correction"
#' Jpn. J. For. Plann. 44:23-39.
#'
#' @param dt_stem_i  Dataframe from the predefined table.
#' @param D  Numeric. DBH in cm.
#' @param H  Numeric. Tree height in m.
#' @param eq_type  Numeric. equation type (currently, 1-4 is defined).
#' @param adj  One of `None`, `3w`, `5w`, `3w_o`, `5w_o`.
#' @return Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#'
#' @importFrom dplyr
#'  %>% filter mutate group_by summarize if_else
#' @export
#'
calcStemVolumeAdj <- function(dt_stem_i, D, H, eq_type = 1, adj = "None"){
  ## switch equations ---
  if(adj == "None"){
    V <- calcStemVolume(dt_stem_i, D, H, eq_type)

  } else if (adj %in% c("3w", "5w")){
    ## DBH is located at the junction of equations?
    th_adj <- if_else(adj == "3w", 3, 5)
    v_bound   <- D - dt_stem_i$D_upper
    wch_close <- abs(v_bound) < th_adj
    is_bound  <- any(wch_close)
    if(!is_bound){
      V <- calcStemVolume(dt_stem_i, D, H, eq_type)
    } else if(adj == "3w"){
      Dt <- dt_stem_i$D_upper[wch_close] - 3
      if(D - Dt < 2){
        V1 <- calcStemVolume(dt_stem_i, Dt,   H, eq_type)
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V2t <- (V1 + V2 + V3)/3
        V <- linearImpute(Dt, V1, Dt+2, V2t, D)
      } else if(D - Dt < 4){
        V1 <- calcStemVolume(dt_stem_i, Dt,   H, eq_type)
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V2t <- (V1 + V2 + V3)/3
        V3t <- (V2 + V3 + V4)/3
        V <- linearImpute(Dt+2, V2t, Dt+4, V3t, D)
      } else {
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V3t <- (V2 + V3 + V4)/3
        V <- linearImpute(Dt+4, V3t, Dt+6, V4, D)
      }

    } else if(adj == "5w"){
      Dt <- dt_stem_i$D_upper[wch_close] - 7
      if(D - Dt < 4){
        V1 <- calcStemVolume(dt_stem_i, Dt,   H, eq_type)
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V5 <- calcStemVolume(dt_stem_i, Dt+8, H, eq_type)
        V3t <- (V1 + V2 + V3 + V4 + V5)/5
        V <- linearImpute(Dt+2, V2, Dt+4, V3t, D)
      } else if(D - Dt < 6){
        V1 <- calcStemVolume(dt_stem_i, Dt,   H, eq_type)
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V5 <- calcStemVolume(dt_stem_i, Dt+8, H, eq_type)
        V6 <- calcStemVolume(dt_stem_i, Dt+10, H, eq_type)
        V3t <- (V1 + V2 + V3 + V4 + V5)/5
        V4t <- (V2 + V3 + V4 + V5 + V6)/5
        V <- linearImpute(Dt+4, V3t, Dt+6, V4t, D)
      } else if(D - Dt < 8){
        V2 <- calcStemVolume(dt_stem_i, Dt+2, H, eq_type)
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V5 <- calcStemVolume(dt_stem_i, Dt+8, H, eq_type)
        V6 <- calcStemVolume(dt_stem_i, Dt+10, H, eq_type)
        V7 <- calcStemVolume(dt_stem_i, Dt+12, H, eq_type)
        V4t <- (V2 + V3 + V4 + V5 + V6)/5
        V5t <- (V3 + V4 + V5 + V6 + V7)/5
        V <- linearImpute(Dt+6, V4t, Dt+8, V5t, D)
      } else if(D - Dt < 10){
        V3 <- calcStemVolume(dt_stem_i, Dt+4, H, eq_type)
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V5 <- calcStemVolume(dt_stem_i, Dt+8, H, eq_type)
        V6 <- calcStemVolume(dt_stem_i, Dt+10, H, eq_type)
        V7 <- calcStemVolume(dt_stem_i, Dt+12, H, eq_type)
        V8 <- calcStemVolume(dt_stem_i, Dt+14, H, eq_type)
        V5t <- (V3 + V4 + V5 + V6 + V7)/5
        V6t <- (V4 + V5 + V6 + V7 + V8)/5
        V <- linearImpute(Dt+8, V5t, Dt+10, V6t, D)
      } else {
        V4 <- calcStemVolume(dt_stem_i, Dt+6, H, eq_type)
        V5 <- calcStemVolume(dt_stem_i, Dt+8, H, eq_type)
        V6 <- calcStemVolume(dt_stem_i, Dt+10, H, eq_type)
        V7 <- calcStemVolume(dt_stem_i, Dt+12, H, eq_type)
        V8 <- calcStemVolume(dt_stem_i, Dt+14, H, eq_type)
        V6t <- (V4 + V5 + V6 + V7 + V8)/5
        V <- linearImpute(Dt+10, V6t, Dt+12, V7, D)
      }
    }

    ### This is not the exact moving average.
  } else if(adj %in% c("3w_o", "5w_o")){
    ## DBH is located at the junction of equations?
    is_bound   <- any(abs(D - dt_stem_i$D_upper) < 1)
    if(!is_bound){
      V <- calcStemVolume(dt_stem_i, D, H, eq_type)
    } else if(adj == "3w_o"){
      V1 <- calcStemVolume(dt_stem_i, D-2, H, eq_type)
      V2 <- calcStemVolume(dt_stem_i, D,   H, eq_type)
      V3 <- calcStemVolume(dt_stem_i, D+2, H, eq_type)
      V <- (V1 + V2 + V3)/3
    } else if(adj == "5w_o"){
      V1 <- calcStemVolume(dt_stem_i, D-4, H, eq_type)
      V2 <- calcStemVolume(dt_stem_i, D-2, H, eq_type)
      V3 <- calcStemVolume(dt_stem_i, D,   H, eq_type)
      V4 <- calcStemVolume(dt_stem_i, D+2, H, eq_type)
      V5 <- calcStemVolume(dt_stem_i, D+4, H, eq_type)
      V <- (V1 + V2 + V3 + V4 + V5)/5
    }

  } else {
    stop("No equation.")
  }

  return(V)
}

#' Calculate stem volume based on junction adjustment, vector version.
#'
#' @description
#' Calculate stem volume using an equation corresponding with the equation type.\cr
#' This function accepts numeric vectors.
#'
#' @details
#' This is a wrapper function to apply 3/5- points moving window average for calculating stem volume when
#' DBH is located at the junction of equations.
#'
#' See the following reference for calculation details:\cr
#' Hosoda et al. (2010):
#' "Differences between the present stem volume tables and the values of the volume equations, and their correction"
#' Jpn. J. For. Plann. 44:23-39.
#'
#' @param dt_stem_i  Dataframe from the predefined table.
#' @param D  Numeric vector. DBH in cm. All elements should be the same DBH class in `dt_stem_i`, but it is acceptable
#'   to pass DBH with different classes (processing is a little slow).
#' @param H  Numeric vector. Tree height in m.
#' @param eq_type  Numeric. equation type (currently, 1-4 is defined).
#' @param adj  Adjusting methods for the junction. One of `None`, `3w`, `5w`, `3w_o`, `5w_o`.
#' @param list_coefs  A list of parameters, which is derived from `dt_stem_i`. Usually, this is set to NULL.
#' @return  Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#'
#' @importFrom dplyr
#'  %>% filter mutate group_by summarize if_else
#' @importFrom purrr
#'  map pmap map_lgl
#' @export
#'
calcStemVolumeMulti <- function(dt_stem_i, D, H, eq_type = 1, adj = "None", list_coefs = NULL){
  ## switch equations ---
  if(adj == "None"){
    V <- calcStemVolume(dt_stem_i, D, H, eq_type, list_coefs)

  } else if (adj %in% c("3w", "5w")){
    ## first calculate tentative V for all observations ---
    V <- calcStemVolume(dt_stem_i, D, H, eq_type, list_coefs)
    ## judge whether DBH is located at the junction of equations.
    ### wch_close: indicates the rows of lower/upper bound (only for ones within lower/upper junctions)
    ### wch_adj: indicates which are bounded within lower/upper junctions
    if(length(D) < length(H) & length(D) == 1) D <- rep(D, length(H)) # adjust if D is shorter
    if(length(D) > length(H) & length(H) == 1) H <- rep(H, length(D)) # adjust if H is shorter
    th_adj <- if_else(adj == "3w", 3, 5)
    row_close <- map(D, function(x) which(abs(x - dt_stem_i$D_upper) < th_adj)) |> unlist()
    wch_adj   <- map_lgl(D, function(x) if_else(any(abs(x - dt_stem_i$D_upper) < th_adj), TRUE, FALSE)) |> which()

    if(adj == "3w"){
      Dt <- dt_stem_i$D_upper[row_close] - 3
      ### condition 1
      is_within <- D[wch_adj] - Dt < 2
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V1 <- calcStemVolume(dt_stem_i, Dt[wch_dt],   H[wch_adj_dt], eq_type, list_coefs)
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V2t <- (V1 + V2 + V3)/3
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt], V1, Dt[wch_dt]+2, V2t, D[wch_adj_dt])
      }
      ### condition 2
      is_within <- D[wch_adj] - Dt < 4 & D[wch_adj] - Dt >= 2
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V1 <- calcStemVolume(dt_stem_i, Dt[wch_dt],   H[wch_adj_dt], eq_type, list_coefs)
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V2t <- (V1 + V2 + V3)/3
        V3t <- (V2 + V3 + V4)/3
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+2, V2t, Dt[wch_dt]+4, V3t, D[wch_adj_dt])
      }
      ### condition 3
      is_within <- D[wch_adj] - Dt >= 4
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V3t <- (V2 + V3 + V4)/3
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+4, V3t, Dt[wch_dt]+6, V4, D[wch_adj_dt])
      }

    } else if(adj == "5w"){
      Dt <- dt_stem_i$D_upper[row_close] - 7
      ### condition 1
      is_within <- D[wch_adj] - Dt < 4
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V1 <- calcStemVolume(dt_stem_i, Dt[wch_dt],   H[wch_adj_dt], eq_type, list_coefs)
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V5 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+8, H[wch_adj_dt], eq_type, list_coefs)
        V3t <- (V1 + V2 + V3 + V4 + V5)/5
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+2, V2, Dt[wch_dt]+4, V3t, D[wch_adj_dt])
      }
      ### condition 2
      is_within <- D[wch_adj] - Dt < 6 & D[wch_adj] - Dt >= 4
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V1 <- calcStemVolume(dt_stem_i, Dt[wch_dt],   H[wch_adj_dt], eq_type, list_coefs)
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V5 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+8, H[wch_adj_dt], eq_type, list_coefs)
        V6 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+10, H[wch_adj_dt], eq_type, list_coefs)
        V3t <- (V1 + V2 + V3 + V4 + V5)/5
        V4t <- (V2 + V3 + V4 + V5 + V6)/5
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+4, V3t, Dt[wch_dt]+6, V4t, D[wch_adj_dt])
      }
      ### condition 3
      is_within <- D[wch_adj] - Dt < 8 & D[wch_adj] - Dt >= 6
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V2 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+2, H[wch_adj_dt], eq_type, list_coefs)
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V5 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+8, H[wch_adj_dt], eq_type, list_coefs)
        V6 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+10, H[wch_adj_dt], eq_type, list_coefs)
        V7 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+12, H[wch_adj_dt], eq_type, list_coefs)
        V4t <- (V2 + V3 + V4 + V5 + V6)/5
        V5t <- (V3 + V4 + V5 + V6 + V7)/5
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+6, V4t, Dt[wch_dt]+8, V5t, D[wch_adj_dt])
      }
      ### condition 4
      is_within <- D[wch_adj] - Dt < 10 & D[wch_adj] - Dt >= 8
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V3 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+4, H[wch_adj_dt], eq_type, list_coefs)
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V5 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+8, H[wch_adj_dt], eq_type, list_coefs)
        V6 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+10, H[wch_adj_dt], eq_type, list_coefs)
        V7 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+12, H[wch_adj_dt], eq_type, list_coefs)
        V8 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+14, H[wch_adj_dt], eq_type, list_coefs)
        V5t <- (V3 + V4 + V5 + V6 + V7)/5
        V6t <- (V4 + V5 + V6 + V7 + V8)/5
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+8, V5t, Dt[wch_dt]+10, V6t, D[wch_adj_dt])
      }
      ### condition 5
      is_within <- D[wch_adj] - Dt >= 10
      if(any(is_within)){
        wch_dt <- which(is_within)
        wch_adj_dt <- wch_adj[wch_dt]
        V4 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+6, H[wch_adj_dt], eq_type, list_coefs)
        V5 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+8, H[wch_adj_dt], eq_type, list_coefs)
        V6 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+10, H[wch_adj_dt], eq_type, list_coefs)
        V7 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+12, H[wch_adj_dt], eq_type, list_coefs)
        V8 <- calcStemVolume(dt_stem_i, Dt[wch_dt]+14, H[wch_adj_dt], eq_type, list_coefs)
        V6t <- (V4 + V5 + V6 + V7 + V8)/5
        V[wch_adj_dt] <- linearImpute(Dt[wch_dt]+10, V6t, Dt[wch_dt]+12, V7, D[wch_adj_dt])
      }
    }

    ### This is not the exact moving average.
  } else if(adj %in% c("3w_o", "5w_o")){
    if(length(D) < length(H) & length(D) == 1) D <- rep(D, length(H)) # adjust if D is shorter
    if(length(D) > length(H) & length(H) == 1) H <- rep(H, length(D)) # adjust if H is shorter
    wch_adj   <- map_lgl(D, function(x) if_else(any(abs(x - dt_stem_i$D_upper) < 1), TRUE, FALSE)) |> which()
     if(adj == "3w_o"){
      V1 <- calcStemVolume(dt_stem_i, D[wch_adj]-2, H[wch_adj], eq_type, list_coefs)
      V2 <- calcStemVolume(dt_stem_i, D[wch_adj],   H[wch_adj], eq_type, list_coefs)
      V3 <- calcStemVolume(dt_stem_i, D[wch_adj]+2, H[wch_adj], eq_type, list_coefs)
      V[wch_adj] <- (V1 + V2 + V3)/3
    }
    if(adj == "5w_o"){
      V1 <- calcStemVolume(dt_stem_i, D[wch_adj]-4, H[wch_adj], eq_type, list_coefs)
      V2 <- calcStemVolume(dt_stem_i, D[wch_adj]-2, H[wch_adj], eq_type, list_coefs)
      V3 <- calcStemVolume(dt_stem_i, D[wch_adj],   H[wch_adj], eq_type, list_coefs)
      V4 <- calcStemVolume(dt_stem_i, D[wch_adj]+2, H[wch_adj], eq_type, list_coefs)
      V5 <- calcStemVolume(dt_stem_i, D[wch_adj]+4, H[wch_adj], eq_type, list_coefs)
      V[wch_adj] <- (V1 + V2 + V3 + V4 + V5)/5
    }

  } else {
    stop("No equation.")
  }

  return(V)
}



#' Calculate stem volume (single)
#'
#' @description
#' Calculate stem volume using coefficients & following the methods described in Hosoda et al. (2010).
#' This function only accepts single numeric value for D & H. To handle numeric vector, use [stemVolume()].\cr
#' VERY SLOW if this is used for processing large amount of data.
#'
#' Different from [stemVolume()], this function can turn off the 3/5 points moving average adjustment if preferred.
#'
#' @details
#' For details, please see the following reference:\cr
#' 細田ら (2010) 現行立木幹材積表と材積式による計算値との相違およびその修正方法. 森林計画学会誌 44: 23-39.
#' (Kazuo HOSODA, Yasushi MITSUDA and Toshiro IEHARA 2010:
#' "Differences between the present stem volume tables and the values of the volume equations, and their correction"
#' Jpn. J. For. Plann. 44:23-39.)
#'
#' There are several differences compared with the description in the original article and
#' computational program developed based on the same article as follows:
#'
#' - 3点移動平均、5点移動平均が複数適用される場合を反映(2024/01/11)
#' - 一部の係数の有効数字で係数の時点でずれが発生している(函館エゾマツ、東京広葉樹の62cm、長野カラマツなど)
#' - 形状高法は係数を直線補間してるが、範囲外も上限/下限で線形補間している。
#' - マクロにおける「札幌トドマツ」、「高知天然スギ」、「青森アカマツ」の修正は反映していない。
#' - 「青森広葉樹」、「高知広葉樹」の5点移動平均の当てはめ範囲がマクロ関数と異なる。
#' - 「北海道針葉樹」の小数第3位への四捨五入で、round()が偶数丸めこみ(五捨五入)を行うため誤差が生じる場合がある。
#'   (2024/01/14)関数を定義して修正したが、内部的な数値のズレで合わない場合が存在する(e.g. H=39mかつD>=81cm)
#'
#' @param Name  Japanese character for identifying the calculation methods (e.g., `東京スギ`).
#'   This is expected from [volumeName()].
#' @param D  Numeric. DBH in cm.
#' @param H  Numeric. Tree height in m.
#' @param list_data  A list of data for coefficients. This is expected from [getStemCoefficients()].
#'   If NULL, internally call [getStemCoefficients()].
#' @param off_adj  Logical. If TRUE, all moving average adjustment for junctions become inactive.
#' @return Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#'
#' @importFrom dplyr %>% filter mutate group_by summarize arrange
#' @importFrom utils tail
#' @export
#'

stemVolumeSingle <- function(Name, D, H, list_data = NULL, off_adj = FALSE){
  ## adjustment to D & H. if either is NA, cause error ---
  if(D < 0) D <- 0
  if(H < 0) H <- 0

  ## read data for coefficients if necessary ---
  if(is.null(list_data)) list_data <- getStemCoefficients()
  dt_stem <- list_data$stem
  dt_ff   <- list_data$ff
  dt_hf   <- list_data$hf

  ## list of Name ---
  ls_stem_name <- c(
    '旭川トドマツ', '旭川エゾマツ', '北見トドマツ', '北見エゾマツ', '帯広トドマツ', '帯広エゾマツ', '札幌トドマツ',
    '札幌エゾマツ', '函館トドマツ', '函館エゾマツ', '函館ヒバ', '函館ブナ', '北海道広葉樹', '青森スギ', '青森アカマツ',
    '青森広葉樹', '秋田スギ', '秋田天然スギ', '秋田カラマツ', '秋田広葉樹', '関東中部モミ', '関東中部ツガ', '前橋スギ',
    '会津新潟スギ', '前橋ヒノキ', '前橋アカマツ', '会津新潟アカマツ', '前橋カラマツ', '前橋広葉樹', '東京スギ',
    '東京ヒノキ', '東京アカマツ', '東京広葉樹', '長野スギ', '長野ヒノキ', '長野天然ヒノキ', '長野サワラ', '長野アカマツ',
    '長野カラマツ', '長野ブナ', '長野サワグルミ', '名古屋スギ', '名古屋ヒノキ', '名古屋天然ヒノキ', '名古屋サワラ',
    '名古屋アカマツ', '名古屋カラマツ', '名古屋広葉樹', '大阪スギ', '山陰天然スギ', '大阪ヒノキ', '大阪アカマツ',
    '山陰アカマツ', '大阪コウヤマキ', '大阪モミ', '大阪広葉樹1型', '大阪広葉樹2型', '高知スギ', '高知天然スギ',
    '高知ヒノキ', '高知天然ヒノキ', '高知アカマツ', '高知モミ', '高知ツガ', '高知広葉樹', '熊本スギ', '飫肥スギ',
    '下屋久天然スギ', '熊本ヒノキ', '熊本モミ', '熊本ツガ', '熊本アカマツ', '霧島アカマツ', '霧島天然アカマツ',
    '熊本広葉樹1類', '熊本広葉樹2類', '上屋久天然スギ'
  )
  ls_FF_name <- c("北海道カラマツ", "北海道針葉樹", "青森天然スギ")
  ls_HF_name <- c("青森ヒバ", "青森針葉樹", "秋田アカマツ")
  ## list of junction adjustment ---
  ls_adjust_3w <- c(
    '旭川トドマツ',  '北見トドマツ', '北見エゾマツ', '帯広トドマツ', '帯広エゾマツ', '札幌トドマツ',
    '函館エゾマツ', '函館ヒバ', '函館ブナ', '北海道広葉樹', '青森アカマツ', '会津新潟アカマツ', '前橋カラマツ',
    '前橋広葉樹', '東京アカマツ', '東京広葉樹', '長野スギ', '名古屋スギ', '名古屋ヒノキ', '名古屋天然ヒノキ',
    '名古屋サワラ', '名古屋広葉樹', '大阪広葉樹1型', '大阪広葉樹2型', '高知スギ', '高知天然スギ',
    '高知ヒノキ', '高知天然ヒノキ', '高知アカマツ', '高知モミ', '高知ツガ', '高知広葉樹', '飫肥スギ',
    '下屋久天然スギ', '上屋久天然スギ', '熊本ヒノキ', '熊本モミ', '熊本ツガ', '熊本アカマツ', '霧島アカマツ')
  ls_adjust_5w <- c('関東中部モミ', '関東中部ツガ', '長野アカマツ', '長野ブナ', '名古屋カラマツ', '大阪アカマツ',
                    '山陰アカマツ', '霧島天然アカマツ', '熊本広葉樹1類', '熊本広葉樹2類')
  ## make inactive adjustment if preferred ---
  if(off_adj){
    ls_adjust_3w <- c("")
    ls_adjust_5w <- c("")
  }

  # 1. stem volume tables -----
  if(Name %in% ls_stem_name){
    ## filter eq. by region+species name
    dt_stem_i <- dt_stem |> filter(name == Name)
    if(nrow(dt_stem_i) == 0) stop("No such Name in equation.")

    ## eq. type 1 -----
    if(dt_stem_i$var_left[1] == "log v" &
       dt_stem_i$var_b[1] == "log d" &
       dt_stem_i$var_c[1] == "log h"){
      # ## get coefs that correspond with DBH size ---
      # l_c <- getCoefs(dt_stem_i, D)
      # ## calculate stem volume ---
      # V <- 10^(l_c[["a"]] +l_c[["b"]]*log10(D) + l_c[["c"]]*log10(H))
      ## calculate stem volume ---
      if(dt_stem_i$name[1] %in% ls_adjust_3w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "3w")
      } else if (dt_stem_i$name[1] %in% ls_adjust_5w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "5w")
      } else {
        V <- calcStemVolume(dt_stem_i, D, H, eq_type = 1)
      }

      ## eq. type 2 -----
    } else if (dt_stem_i$var_left[1] == "v" &
               dt_stem_i$var_a[1] == "h" &
               dt_stem_i$var_c[1] == "d2h"){
      # ## get coefs that correspond with DBH size ---
      # l_c <- getCoefs(dt_stem_i, D)
      # ## DBH unit conversion if required
      # if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      # ## calculate stem volume ---
      # l_c$b[is.na(l_c$b)] <- 0 # for Kochi akamatsu
      # V <- l_c[["a"]]*H + l_c[["b"]]*Dm^2 + l_c[["c"]]*Dm^2*H + l_c[["d"]]
      ## calculate stem volume ---
      if(dt_stem_i$name[1] %in% ls_adjust_3w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 2, adj = "3w")
      } else if (dt_stem_i$name[1] %in% ls_adjust_5w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 2, adj = "5w")
      } else {
        V <- calcStemVolume(dt_stem_i, D, H, eq_type = 2)
      }

      ## eq. type 3 -----
    } else if (dt_stem_i$var_left[1] == "v" &
               dt_stem_i$var_b[1] == "" &
               dt_stem_i$var_c[1] == "d2h"){
      # ## get coefs that correspond with DBH size ---
      # l_c <- getCoefs(dt_stem_i, D)
      # ## DBH unit conversion if required
      # if(dt_stem_i$D_unit[1] == "m") Dm <- D / 100 # cm to m
      # ## calculate stem volume ---
      # V <- l_c[["c"]]*Dm^2*H + l_c[["d"]]
      ## calculate stem volume ---
      if(dt_stem_i$name[1] %in% ls_adjust_3w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 3, adj = "3w")
      } else if (dt_stem_i$name[1] %in% ls_adjust_5w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 3, adj = "5w")
      } else {
        V <- calcStemVolume(dt_stem_i, D, H, eq_type = 3)
      }

      ## eq. type 4 -----
    } else if (dt_stem_i$var_left[1] == "log v" &
               dt_stem_i$var_b[1] == "d2h" &
               dt_stem_i$var_c[1] == ""){
      # ## get coefs that correspond with DBH size ---
      # l_c <- getCoefs(dt_stem_i, D)
      # ## calculate stem volume ---
      # V <- 10^(l_c[["a"]] + l_c[["b"]]*log10(D^2*H))
      ## calculate stem volume ---
      if(dt_stem_i$name[1] %in% ls_adjust_3w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 4, adj = "3w")
      } else if (dt_stem_i$name[1] %in% ls_adjust_5w){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 4, adj = "5w")
      } else {
        V <- calcStemVolume(dt_stem_i, D, H, eq_type = 4)
      }

    } else {
      stop("No such equation.")
    }

    ## special adjustment -----
    if(Name == "上屋久天然スギ" & D >= 61){
      V = 0.1513 + 0.8898 * V
    }
    if(!off_adj){
      if(Name == "上屋久天然スギ" & abs(D - 61) < 5){
        dt_stem_s <- dt_stem_i |> rbind(dt_stem_i[2,]) |> arrange(D_lower)
        dt_stem_s$D_upper[2] <- 61
        dt_stem_s$D_lower[3] <- 61
        V <- calcStemVolumeAdj(dt_stem_s, D, H, eq_type = 5, adj = "5w")
      }
      if(Name == "札幌トドマツ" & D >= 95){
        dt_stem_s <- dt_stem |> filter(name == Name) |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        V <- calcStemVolumeAdj(dt_stem_s, D, H, eq_type = 1, adj = "5w")
      }
      if(Name == "札幌トドマツ" & abs(D - 95) < 5){
        dt_stem_s <- dt_stem |> filter(name == Name) |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        V <- calcStemVolumeAdj(dt_stem_s, D, H, eq_type = 1, adj = "5w")
      }
      if(Name == "青森アカマツ" & abs(D - 47) < 5){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "5w")
      }
      if(Name == "青森広葉樹" & abs(D - 71) < 5){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "5w")
      }
      if(Name == "高知天然スギ" & abs(D - 101) < 5){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "5w")
      }
      if(Name == "高知広葉樹" & abs(D - 61) < 5){
        V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 2, adj = "5w")
      }
      if(Name == "東京スギ" & H < 11.5 & D >= 36 & D <= 70){
        V <- linearImputeVolume(dt_stem_i, 36, 70, D, H, eq_type = 1)
      } else if(Name == "東京スギ" & H < 12.5 & D >= 36 & D <= 68){
        V <- linearImputeVolume(dt_stem_i, 36, 68, D, H, eq_type = 1)
      } else if(Name == "東京スギ" & H < 14.5 & D >= 36 & D <= 62){
        V <- linearImputeVolume(dt_stem_i, 36, 62, D, H, eq_type = 1)
      } else if(Name == "東京スギ" & H < 16.5 & D >= 36 & D <= 60){
        V <- linearImputeVolume(dt_stem_i, 36, 60, D, H, eq_type = 1)
      } else if(Name == "東京スギ" & H < 18.5 & D >= 40 & D <= 60){
        v1 <- calcStemVolumeAdj(dt_stem_i, 40, H, eq_type = 1, adj = "3w")
        v2 <- calcStemVolume(dt_stem_i, 60, H, eq_type = 1)
        V  <- linearImpute(40, v1, 60, v2, D)
      } else if(Name == "東京スギ" & H < 19.5 & D >= 40 & D <= 58){
        v1 <- calcStemVolumeAdj(dt_stem_i, 40, H, eq_type = 1, adj = "3w")
        v2 <- calcStemVolume(dt_stem_i, 58, H, eq_type = 1)
        V  <- linearImpute(40, v1, 58, v2, D)
      } else if(Name == "東京スギ" & H < 20.5 & D >= 40 & D <= 54){
        v1 <- calcStemVolumeAdj(dt_stem_i, 40, H, eq_type = 1, adj = "3w")
        v2 <- calcStemVolume(dt_stem_i, 54, H, eq_type = 1)
        V  <- linearImpute(40, v1, 54, v2, D)
      } else if(Name == "東京スギ" & H < 21.5 & D >= 40 & D <= 52){
        v1 <- calcStemVolumeAdj(dt_stem_i, 40, H, eq_type = 1, adj = "3w")
        v2 <- calcStemVolume(dt_stem_i, 52, H, eq_type = 1)
        V  <- linearImpute(40, v1, 52, v2, D)
      } else if(Name == "東京スギ" & H < 22.5 & D >= 40 & D <= 50){
        v1 <- calcStemVolumeAdj(dt_stem_i, 40, H, eq_type = 1, adj = "3w")
        v2 <- calcStemVolume(dt_stem_i, 50, H, eq_type = 1)
        V  <- linearImpute(40, v1, 50, v2, D)
      } else {
        V <- V
      }
    }



  # 2. DBH form factor -----
  } else if (Name %in% ls_FF_name){
    ## get f coefficient ---
    if(Name %in% unique(dt_ff$name)){
      D_f <- floor(D/2)*2
      f1 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == D_f]
      f2 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == D_f + 2]
      ## adjust if out of range ---
      if(D >= max(dt_ff$DBH)){
        f1 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == max(dt_ff$DBH) - 2]
        f2 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == max(dt_ff$DBH)]
        D_f <- max(dt_ff$DBH) - 2
      }
      if(D < min(dt_ff$DBH)){
        f1 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == min(dt_ff$DBH)]
        f2 <- dt_ff$factor[dt_ff$name == Name & dt_ff$DBH == min(dt_ff$DBH) + 2]
        D_f <- min(dt_ff$DBH)
      }
      f <- linearImpute(D_f, f1, D_f+2, f2, D)
    } else if(Name == "北海道カラマツ") {
      f_d <- 0.439004 + 0.916461/D - 0.073809/D^2
      f_h <- 0.435719 + 0.515867/H + 2.481278/H^2
      f <- (round2(f_d, 3) + round2(f_h, 3))/2
    } else if(Name == "北海道針葉樹") {
      f_d <- 0.5 - 0.0008*D + 0.421*exp(-0.12*D)
      f_h <- 0.61 - 0.0055*H + 5.48*exp(-1.025*H)
      f <- (round2(f_d, 3) + round2(f_h, 3))/2
    }
    ## calculate stem volume ---
    V <- H * D^2 * pi * f / 40000
    V[is.na(V)] <- 0

  # 3. Height form -----
  } else if (Name %in% ls_HF_name){
    ## get hf coefficient ---
    dt_hf_i <- dt_hf |> filter(name == Name)
    H_hf_i  <- floor(H)
    hf1 <- dt_hf_i$factor[dt_hf_i$H == H_hf_i]
    hf2 <- dt_hf_i$factor[dt_hf_i$H == H_hf_i+1]
    ## adjust if out of range ---
    if(H >= max(dt_hf_i$H)){
      hf1 <- dt_hf_i$factor[dt_hf_i$H == max(dt_hf_i$H) - 1]
      hf2 <- dt_hf_i$factor[dt_hf_i$H == max(dt_hf_i$H)]
      H_hf_i <- max(dt_hf_i$H) - 1
    }
    if(H < min(dt_hf_i$H)){
      hf1 <- dt_hf_i$factor[dt_hf_i$H == min(dt_hf_i$H)]
      hf2 <- dt_hf_i$factor[dt_hf_i$H == min(dt_hf_i$H) + 1]
      H_hf_i <-  min(dt_hf_i$H)
    }
    ## calculate stem volume ---
    hf <- linearImpute(H_hf_i, hf1, H_hf_i+1, hf2, H)
    V <- D^2 * pi * hf / 40000

  } else {
    stop("No such Name in equation.")
  }

  return(V)
}



#' Calculate stem volume
#'
#' @description
#' Calculate stem volume using coefficients & following the methods described in Hosoda et al. (2010).
#' This function accepts numeric vectors for D & H.
#'
#' @details
#' For details, please see the following reference:\cr
#' 細田ら (2010) 現行立木幹材積表と材積式による計算値との相違およびその修正方法. 森林計画学会誌 44: 23-39.\cr
#' (Kazuo HOSODA, Yasushi MITSUDA and Toshiro IEHARA 2010:
#' "Differences between the present stem volume tables and the values of the volume equations, and their correction"
#' Jpn. J. For. Plann. 44:23-39.)
#'
#' There are several differences compared with the description in the original article and
#' computational program developed based on the same article as follows:
#'
#' - 3点移動平均、5点移動平均が複数適用される場合を反映(2024/01/11)
#' - 一部の係数の有効数字で係数の時点でずれが発生している(函館エゾマツ、東京広葉樹の62cm、長野カラマツなど)
#' - 形状高法は係数を直線補間してるが、範囲外も上限/下限で線形補間している。
#' - マクロにおける「札幌トドマツ」、「高知天然スギ」、「青森アカマツ」の修正は反映していない。
#' - 「青森広葉樹」、「高知広葉樹」の5点移動平均の当てはめ範囲がマクロ関数と異なる。
#' - 「北海道針葉樹」の小数第3位への四捨五入で、round()が偶数丸めこみ(五捨五入)を行うため誤差が生じる場合がある。
#'   (2024/01/14)関数を定義して修正したが、内部的な数値のズレで合わない場合が存在する(e.g. H=39mかつD>=81cm)
#'
#' @param Name  Japanese character for identifying the calculation methods (e.g., `東京スギ`).
#'   This is expected from [volumeName()].
#' @param D  Numeric. DBH in cm.
#' @param H  Numeric. Tree height in m.
#' @return  Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#' @importFrom dplyr %>% filter mutate group_by summarize arrange
#' @importFrom purrr map
#' @importFrom utils tail
#' @export
#'
stemVolume <- function(Name, D, H){
  ## check the length ---
  if(length(Name) == 1 & length(D) > 1) Name <- rep(Name, length(D)) # in case single Name is provided.
  if(!all(length(Name) == c(length(D), length(H)))) stop("Different data length.")
  ## adjustment to D & H. if either is NA, cause error ---
  if(any(D < 0)) D[D < 0] <- 0
  if(any(H < 0)) H[H < 0] <- 0
  V <- rep(NA, length(D))

  ## data for coefficients ---
  list_data <- getStemCoefficients()
  dt_stem <- list_data$stem
  dt_ff   <- list_data$ff
  dt_hf   <- list_data$hf


  ## list of Name ---
  ls_stem_name <- c(
    '旭川トドマツ', '旭川エゾマツ', '北見トドマツ', '北見エゾマツ', '帯広トドマツ', '帯広エゾマツ', '札幌トドマツ',
    '札幌エゾマツ', '函館トドマツ', '函館エゾマツ', '函館ヒバ', '函館ブナ', '北海道広葉樹', '青森スギ', '青森アカマツ',
    '青森広葉樹', '秋田スギ', '秋田天然スギ', '秋田カラマツ', '秋田広葉樹', '関東中部モミ', '関東中部ツガ', '前橋スギ',
    '会津新潟スギ', '前橋ヒノキ', '前橋アカマツ', '会津新潟アカマツ', '前橋カラマツ', '前橋広葉樹', '東京スギ',
    '東京ヒノキ', '東京アカマツ', '東京広葉樹', '長野スギ', '長野ヒノキ', '長野天然ヒノキ', '長野サワラ', '長野アカマツ',
    '長野カラマツ', '長野ブナ', '長野サワグルミ', '名古屋スギ', '名古屋ヒノキ', '名古屋天然ヒノキ', '名古屋サワラ',
    '名古屋アカマツ', '名古屋カラマツ', '名古屋広葉樹', '大阪スギ', '山陰天然スギ', '大阪ヒノキ', '大阪アカマツ',
    '山陰アカマツ', '大阪コウヤマキ', '大阪モミ', '大阪広葉樹1型', '大阪広葉樹2型', '高知スギ', '高知天然スギ',
    '高知ヒノキ', '高知天然ヒノキ', '高知アカマツ', '高知モミ', '高知ツガ', '高知広葉樹', '熊本スギ', '飫肥スギ',
    '下屋久天然スギ', '熊本ヒノキ', '熊本モミ', '熊本ツガ', '熊本アカマツ', '霧島アカマツ', '霧島天然アカマツ',
    '熊本広葉樹1類', '熊本広葉樹2類', '上屋久天然スギ'
  )
  ls_FF_name <- c("北海道カラマツ", "北海道針葉樹", "青森天然スギ")
  ls_HF_name <- c("青森ヒバ", "青森針葉樹", "秋田アカマツ")
  ## list of junction adjustment ---
  ls_adjust_3w <- c(
    '旭川トドマツ',  '北見トドマツ', '北見エゾマツ', '帯広トドマツ', '帯広エゾマツ', '札幌トドマツ',
    '函館エゾマツ', '函館ヒバ', '函館ブナ', '北海道広葉樹', '青森アカマツ', '会津新潟アカマツ', '前橋カラマツ',
    '前橋広葉樹', '東京アカマツ', '東京広葉樹', '長野スギ', '名古屋スギ', '名古屋ヒノキ', '名古屋天然ヒノキ',
    '名古屋サワラ', '名古屋広葉樹', '大阪広葉樹1型', '大阪広葉樹2型', '高知スギ', '高知天然スギ',
    '高知ヒノキ', '高知天然ヒノキ', '高知アカマツ', '高知モミ', '高知ツガ', '高知広葉樹', '飫肥スギ',
    '下屋久天然スギ', '上屋久天然スギ', '熊本ヒノキ', '熊本モミ', '熊本ツガ', '熊本アカマツ', '霧島アカマツ')
  ls_adjust_5w <- c('関東中部モミ', '関東中部ツガ', '長野アカマツ', '長野ブナ', '名古屋カラマツ', '大阪アカマツ',
                    '山陰アカマツ', '霧島天然アカマツ', '熊本広葉樹1類', '熊本広葉樹2類')


  ## loop over each unique Name ---
  ls_Name <- unique(Name)
  if(!all(ls_Name %in% c(ls_stem_name, ls_FF_name, ls_HF_name)))
    stop("There are Names that are not listed in this caclulation.")

  for(i in 1:length(ls_Name)){
    ## extract D & H ---
    wch_i <- which(Name == ls_Name[i])
    D_i <- D[wch_i]
    H_i <- H[wch_i]
    V_i <- rep(NA, length(wch_i))

    # 1. stem volume tables -----
    if(ls_Name[i] %in% ls_stem_name){
      ## get coefficients ---
      dt_stem_i <- dt_stem |> filter(name == ls_Name[i])
      ## get equation type ---
      if(dt_stem_i$var_left[1] == "log v" & dt_stem_i$var_b[1] == "log d" & dt_stem_i$var_c[1] == "log h"){
        eq_type_i <- 1
      } else if (dt_stem_i$var_left[1] == "v" & dt_stem_i$var_a[1] == "h" & dt_stem_i$var_c[1] == "d2h"){
        eq_type_i <- 2
      } else if (dt_stem_i$var_left[1] == "v" & dt_stem_i$var_b[1] == "" & dt_stem_i$var_c[1] == "d2h"){
        eq_type_i <- 3
      } else if (dt_stem_i$var_left[1] == "log v" & dt_stem_i$var_b[1] == "d2h" & dt_stem_i$var_c[1] == ""){
        eq_type_i <- 4
      } else {
        stop("No such equation.")
      }

      ## loop for each DBH class (with/without adjustment) ---
      ## NOTE: This can be done without for-loop because calcStemVolumeMulti() and calcStemVolume()
      ## can handle different DBH classes, but it is slower than a process with for-loop.
      for(j in 1:nrow(dt_stem_i)){
        wch_j <- which(dt_stem_i$D_lower[j] <= D_i & dt_stem_i$D_upper[j] > D_i)
        if(length(wch_j) == 0) next
        D_j <- D_i[wch_j]
        H_j <- H_i[wch_j]

        ## calculate stem volume ---
        ### switch for adjustment ---
        if(ls_Name[i] %in% ls_adjust_3w){
          ### 3 points moving windows
          V_j <- calcStemVolumeMulti(dt_stem_i, D_j, H_j, eq_type_i, adj = "3w")
        } else if (ls_Name[i] %in% ls_adjust_5w){
          ### 5 points moving windows
          V_j <- calcStemVolumeMulti(dt_stem_i, D_j, H_j, eq_type_i, adj = "5w")
        } else {
          ### no adjustment
          V_j <- calcStemVolume(dt_stem_i, D_j, H_j, eq_type_i)
        }
        ## combine Volume for loop-i ---
        V_i[wch_j] <- V_j
      }

      ## special adjustment -----
      if(ls_Name[i] == "上屋久天然スギ" & any(D_i >= 61)){
        wch_replace <- which(D_i >= 61)
        V_i[wch_replace] = 0.1513 + 0.8898 * V_i[wch_replace]
      }
      if(ls_Name[i] == "上屋久天然スギ" & any(abs(D_i - 61) < 5)){
        dt_stem_s <- dt_stem_i |> rbind(dt_stem_i[2,]) |> arrange(D_lower)
        dt_stem_s$D_upper[2] <- 61
        dt_stem_s$D_lower[3] <- 61
        wch_replace <- which(abs(D_i - 61) < 5)
        V_i[wch_replace] = calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 5, adj = "5w")
      }
      if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 95)){
        dt_stem_s <- dt_stem |> filter(name == "札幌トドマツ") |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        wch_replace <- which(D_i >= 95)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "札幌トドマツ" & any(abs(D_i - 95) < 5)){
        dt_stem_s <- dt_stem |> filter(name == "札幌トドマツ") |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        wch_replace <- which(abs(D_i - 95) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "青森アカマツ" & any(abs(D_i - 47) < 5)){
        wch_replace <- which(abs(D_i - 47) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "青森広葉樹" & any(abs(D_i - 71) < 5)){
        wch_replace <- which(abs(D_i - 71) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "高知天然スギ" & any(abs(D_i - 101) < 5)){
        wch_replace <- which(abs(D_i - 101) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "高知広葉樹" & any(abs(D_i - 61) < 5)){
        wch_replace <- which(abs(D_i - 61) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 2, adj = "5w")
      }
      if(ls_Name[i] == "東京スギ" & any(H_i < 11.5 & D_i >= 36 & D_i <= 70)){
        wch_replace <- which(H_i < 11.5 & D_i >= 36 & D_i <= 70)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 70, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 11.5 & H_i < 12.5 & D_i >= 36 & D_i <= 68)){
        wch_replace <- which(H_i >= 11.5 & H_i < 12.5 & D_i >= 36 & D_i <= 68)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 68, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 12.5 & H_i < 14.5 & D_i >= 36 & D_i <= 62)){
        wch_replace <- which(H_i >= 12.5 & H_i < 14.5 & D_i >= 36 & D_i <= 62)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 62, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 14.5 & H_i < 16.5 & D_i >= 36 & D_i <= 60)){
        wch_replace <- which(H_i >= 14.5 & H_i < 16.5 & D_i >= 36 & D_i <= 60)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 60, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 16.5 & H_i < 18.5 & D_i >= 40 & D_i <= 60)){
        wch_replace <- which(H_i >= 16.5 & H_i < 18.5 & D_i >= 40 & D_i <= 60)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 60, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 60, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 18.5 & H_i < 19.5 & D_i >= 40 & D_i <= 58)){
        wch_replace <- which(H_i >= 18.5 & H_i < 19.5 & D_i >= 40 & D_i <= 58)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 58, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 58, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 19.5 & H_i < 20.5 & D_i >= 40 & D_i <= 54)){
        wch_replace <- which(H_i >= 19.5 & H_i < 20.5 & D_i >= 40 & D_i <= 54)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 54, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 54, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 20.5 & H_i < 21.5 & D_i >= 40 & D_i <= 52)){
        wch_replace <- which(H_i >= 20.5 & H_i < 21.5 & D_i >= 40 & D_i <= 52)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 52, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 52, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 21.5 & H_i < 22.5 & D_i >= 40 & D_i <= 50)){
        wch_replace <- which(H_i >= 21.5 & H_i < 22.5 & D_i >= 40 & D_i <= 50)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 50, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 50, v2, D_i[wch_replace])
      }


      # 2. DBH form factor -----
    } else if (ls_Name[i] %in% ls_FF_name){
      ## get f coefficient ---
      if(ls_Name[i] %in% unique(dt_ff$name)){
        D_f <- floor(D_i/2)*2
        f1 <- map(D_i, function(x) dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == floor(x/2) * 2])
        f1[lengths(f1) == 0] <- NA
        f1 <- unlist(f1)
        f2 <- map(D_i, function(x) dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == floor(x/2) * 2 + 2])
        f2[lengths(f2) == 0] <- NA
        f2 <- unlist(f2)
        ## adjust if out of range ---
        if(any(D_i >= max(dt_ff$DBH))){
          wch_replace <- which(D_i >= max(dt_ff$DBH))
          f1[wch_replace] <- dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == max(dt_ff$DBH) - 2]
          f2[wch_replace] <- dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == max(dt_ff$DBH)]
          D_f[wch_replace] <- max(dt_ff$DBH) - 2
        }
        if(any(D_i < min(dt_ff$DBH))){
          wch_replace <- which(D_i < min(dt_ff$DBH))
          f1[wch_replace] <- dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == min(dt_ff$DBH)]
          f2[wch_replace] <- dt_ff$factor[dt_ff$name == ls_Name[i] & dt_ff$DBH == min(dt_ff$DBH) + 2]
          D_f[wch_replace] <- min(dt_ff$DBH)
        }
        f <- linearImpute(D_f, f1, D_f+2, f2, D_i)
      } else if(ls_Name[i] == "北海道カラマツ") {
        f_d <- 0.439004 + 0.916461/D_i- 0.073809/D_i^2
        f_h <- 0.435719 + 0.515867/H_i + 2.481278/H_i^2
        f <- (round2(f_d, 3) + round2(f_h, 3))/2
      } else if(ls_Name[i] == "北海道針葉樹") {
        f_d <- 0.5 - 0.0008*D_i+ 0.421*exp(-0.12*D_i)
        f_h <- 0.61 - 0.0055*H_i + 5.48*exp(-1.025*H_i)
        f <- (round2(f_d, 3) + round2(f_h, 3))/2
      }
      ## calculate stem volume ---
      V_i <- H_i * D_i^2 * pi * f / 40000
      V_i[is.na(V_i)] <- 0  # in case D = 0 or H = 0

      # 3. Height form -----
    } else if (ls_Name[i] %in% ls_HF_name){
      ## get hf coefficient ---
      dt_hf_i <- dt_hf |> filter(name == ls_Name[i])
      H_hf_i  <- floor(H_i)
      hf1 <- map(H_i, function(x) dt_hf_i$factor[dt_hf_i$H == floor(x)])
      hf1[lengths(hf1) == 0] <- NA
      hf1 <- unlist(hf1)
      hf2 <- map(H_i, function(x) dt_hf_i$factor[dt_hf_i$H == floor(x)+1])
      hf2[lengths(hf2) == 0] <- NA
      hf2 <- unlist(hf2)
      ## adjust if out of range ---
      if(any(H_i >= max(dt_hf_i$H))){
        wch_replace <- which(H_i >= max(dt_hf_i$H))
        hf1[wch_replace] <- dt_hf_i$factor[dt_hf_i$H == max(dt_hf_i$H) - 1]
        hf2[wch_replace] <- dt_hf_i$factor[dt_hf_i$H == max(dt_hf_i$H)]
        H_hf_i[wch_replace] <- max(dt_hf_i$H) - 1
      }
      if(any(H_i < min(dt_hf_i$H))){
        wch_replace <- which(H_i < min(dt_hf_i$H))
        hf1[wch_replace] <- dt_hf_i$factor[dt_hf_i$H == min(dt_hf_i$H)]
        hf2[wch_replace] <- dt_hf_i$factor[dt_hf_i$H == min(dt_hf_i$H) + 1]
        H_hf_i[wch_replace] <- min(dt_hf_i$H)
      }
      ## calculate stem volume ---
      hf <- linearImpute(H_hf_i, hf1, H_hf_i+1, hf2, H_i)
      V_i <- D_i^2 * pi * hf / 40000

    } else {
      stop("No such Name in equation.")
    }

    ## combine the calculated volume ---
    V[wch_i] <- V_i
  }

  return(V)
}




#' Get Region name (single)
#'
#' Get Region name that is used in [stemVolume()].
#' Only accept single character for Region & Spp. To handle vector, use a wrapper function [volumeName()].
#'
#' Species & Region names are defined following Hosoda et al. (2010) Jpn. J. For. Plann. 44:23-39.
#'  (https://doi.org/10.20659/jjfp.44.2_23) and
#'  "幹材積計算プログラム" developed by FFPRI.
#' (https://www.ffpri.affrc.go.jp/database/stemvolume/index.html)
#'
#' @param Region  Character. Region/Prefecture name that are listed in `RS`.
#'   Please see the original article by Hosoda et al. (2010) and description in the above program.
#' @param Spp  Character. Species name in Japanese. If it does not match any species in the lists, then `広葉樹` is returned.
#'   Please see the original article by Hosoda et al. (2010) and description in the above program.
#' @param RS  A list of parameters, which is generated by [getRegionName()].
#'   List names should be strictly the same as one generated by [getRegionName()] to derive the correct parameters.
#' @param name_invalid  Character to determine the behavior when invalid region name is provided for `Region`.
#'   If NULL, cause error with message. If character string provided, return that string for invalid ones.
#' @return  Name for identifying and calculating stem volume in [stemVolume()]
#'
#' @importFrom dplyr %>% filter mutate group_by summarize
#' @importFrom stringr str_replace str_detect
#' @importFrom stringi stri_trans_general
#' @export

volumeNameSingle <- function(Region, Spp, RS = NULL, name_invalid = NULL){
  # prepare for Region ---
  name_region <- Region
  name_region <- str_replace(name_region, "(府|県)$", "")
  name_region <- str_replace(name_region, "東京都", "東京") # to avoid "京都" to "都"
  # prepare for Spp ---
  name_spp <- Spp
  if(str_detect(name_spp, pattern = "\\p{Katakana}|\\p{Han}"))
    name_spp <- stri_trans_general(name_spp, "Halfwidth-Fullwidth")
  # prepare for list RS ---
  if(is.null(RS)) RS <- getRegionName()


  ## 1. Asahikawa
  if(name_region %in% RS$r_1_Asahikawa){
    if(name_spp %in% RS$s_1_Todomatsu){
      name_RS <- "旭川トドマツ"
    } else if(name_spp %in% RS$s_1_Ezomatsu){
      name_RS <- "旭川エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% RS$s_x_Sugi){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 2. Kitami
  } else if (name_region %in% RS$r_2_Kitami){
    if(name_spp %in% RS$s_2_Todomatsu){
      name_RS <- "北見トドマツ"
    } else if(name_spp %in% RS$s_2_Ezomatsu){
      name_RS <- "北見エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% RS$s_x_Sugi){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 3. Obihiro
  } else if (name_region %in% RS$r_3_Obihiro){
    if(name_spp %in% RS$s_3_Todomatsu){
      name_RS <- "帯広トドマツ"
    } else if(name_spp %in% RS$s_3_Ezomatsu){
      name_RS <- "帯広エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% RS$s_x_Sugi){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 4. Sapporo
  } else if (name_region %in% RS$r_4_Sapporo){
    if(name_spp %in% RS$s_4_Todomatsu){
      name_RS <- "札幌トドマツ"
    } else if(name_spp %in% RS$s_4_Ezomatsu){
      name_RS <- "札幌エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% RS$s_x_Sugi){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 5. Hakodate
  } else if (name_region %in% RS$r_5_Hakodate){
    if(name_spp %in% RS$s_5_Todomatsu){
      name_RS <- "函館トドマツ"
    } else if(name_spp %in% RS$s_5_Ezomatsu){
      name_RS <- "函館エゾマツ"
    } else if(name_spp %in% RS$s_5_Hiba){
      name_RS <- "函館ヒバ"
    } else if(name_spp %in% RS$s_5_Buna){
      name_RS <- "函館ブナ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% RS$s_x_Sugi){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 6. Aomori
  } else if (name_region %in% RS$r_6_Aomori){
    if(name_spp %in% RS$s_6_Sugi){
      name_RS <- "青森スギ"
    } else if(name_spp %in% RS$s_6_TenSugi){
      name_RS <- "青森天然スギ"
    } else if(name_spp %in% RS$s_6_Hiba){
      name_RS <- "青森ヒバ"
    } else if(name_spp %in% RS$s_6_Akamatsu){
      name_RS <- "青森アカマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "秋田カラマツ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "青森針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "青森広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 7. Akita
  } else if (name_region %in% RS$r_7_Akita){
    if(name_spp %in% RS$s_7_Sugi){
      name_RS <- "秋田スギ"
    } else if(name_spp %in% RS$s_7_TenSugi){
      name_RS <- "秋田天然スギ"
    } else if(name_spp %in% RS$s_7_Karamatsu){
      name_RS <- "秋田カラマツ"
    } else if(name_spp %in% RS$s_7_Akamatsu){
      name_RS <- "秋田アカマツ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "青森針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "秋田広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 8. Maebashi
  } else if (name_region %in% RS$r_8_Maebashi){
    if(name_spp %in% RS$s_8_Sugi & name_region %in% c("会津新潟", "新潟")){
      name_RS <- "会津新潟スギ"
    } else if(name_spp %in% RS$s_8_Sugi){
      name_RS <- "前橋スギ"
    } else if(name_spp %in% RS$s_8_Hinoki){
      name_RS <- "前橋ヒノキ"
    } else if(name_spp %in% RS$s_8_Karamatsu){
      name_RS <- "前橋カラマツ"
    } else if(name_spp %in% RS$s_8_Akamatsu & name_region %in% c("会津新潟", "新潟")){
      name_RS <- "会津新潟アカマツ"
    } else if(name_spp %in% RS$s_8_Akamatsu){
      name_RS <- "前橋アカマツ"
    } else if(name_spp %in% RS$s_k_Tsuga){
      name_RS <- "関東中部ツガ"
    } else if(name_spp %in% RS$s_k_Momi){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "前橋広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 9. Tokyo
  } else if (name_region %in% RS$r_9_Tokyo){
    if(name_spp %in% RS$s_9_Sugi){
      name_RS <- "東京スギ"
    } else if(name_spp %in% RS$s_9_Hinoki){
      name_RS <- "東京ヒノキ"
    } else if(name_spp %in% RS$s_9_Akamatsu){
      name_RS <- "東京アカマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "長野カラマツ"
    } else if(name_spp %in% RS$s_k_Tsuga){
      name_RS <- "関東中部ツガ"
    } else if(name_spp %in% RS$s_k_Momi){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "東京広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 10. Nagano
  } else if (name_region %in% RS$r_10_Nagano){
    if(name_spp %in% RS$s_10_Sugi){
      name_RS <- "長野スギ"
    } else if(name_spp %in% RS$s_10_Hinoki){
      name_RS <- "長野ヒノキ"
    } else if(name_spp %in% RS$s_10_TenHinoki){
      name_RS <- "長野天然ヒノキ"
    } else if(name_spp %in% RS$s_10_Sawara){
      name_RS <- "長野サワラ"
    } else if(name_spp %in% RS$s_10_Akamatsu){
      name_RS <- "長野アカマツ"
    } else if(name_spp %in% RS$s_10_Karamatsu){
      name_RS <- "長野カラマツ"
    } else if(name_spp %in% RS$s_k_Tsuga){
      name_RS <- "関東中部ツガ"
    } else if(name_spp %in% RS$s_k_Momi){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "関東中部モミ"
    } else if(name_spp %in% RS$s_10_Sawagurumi){
      name_RS <- "長野サワグルミ"
    } else if(name_spp %in% RS$s_10_Buna){
      name_RS <- "長野ブナ"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "長野ブナ"
    } else {
      stop("No such species in this region.")
    }
    ## 11. Nagoya
  } else if (name_region %in% RS$r_11_Nagoya){
    if(name_spp %in% RS$s_11_Sugi){
      name_RS <- "名古屋スギ"
    } else if(name_spp %in% RS$s_11_Hinoki){
      name_RS <- "名古屋ヒノキ"
    } else if(name_spp %in% RS$s_11_TenHinoki){
      name_RS <- "名古屋天然ヒノキ"
    } else if(name_spp %in% RS$s_11_Sawara){
      name_RS <- "名古屋サワラ"
    } else if(name_spp %in% RS$s_11_Akamatsu){
      name_RS <- "名古屋アカマツ"
    } else if(name_spp %in% RS$s_11_Karamatsu){
      name_RS <- "名古屋カラマツ"
    } else if(name_spp %in% RS$s_k_Tsuga){
      name_RS <- "関東中部ツガ"
    } else if(name_spp %in% RS$s_k_Momi){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "関東中部モミ"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "名古屋広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 12. Osaka
  } else if (name_region %in% RS$r_12_Osaka){
    if(name_spp %in% RS$s_12_Sugi){
      name_RS <- "大阪スギ"
    } else if(name_spp %in% RS$s_12_TenSugi){
      name_RS <- "山陰天然スギ"
    } else if(name_spp %in% RS$s_12_Hinoki){
      name_RS <- "大阪ヒノキ"
    } else if(name_spp %in% RS$s_12_Akamatsu & name_region %in% c("山陰", "鳥取", "島根", "石川", "福井",
                                                                  "中丹地方", "丹後地方", "但馬地方")){
      name_RS <- "山陰アカマツ"
    } else if(name_spp %in% RS$s_12_Akamatsu){
      name_RS <- "大阪アカマツ"
    } else if(name_spp %in% RS$s_12_Kouyamaki){
      name_RS <- "大阪コウヤマキ"
    } else if(name_spp %in% RS$s_12_Momi){
      name_RS <- "大阪モミ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "名古屋カラマツ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "山陰アカマツ"
    } else if(any(str_detect(name_spp, RS$s_12_Broad_1))){
      name_RS <- "大阪広葉樹1型"
    } else if(any(str_detect(name_spp, RS$s_12_Broad_2))){
      name_RS <- "大阪広葉樹2型"
    } else {
      stop("No such species in this region.")
    }
    ## 13. Kochi
  } else if (name_region %in% RS$r_13_Kochi){
    if(name_spp %in% RS$s_13_Sugi){
      name_RS <- "高知スギ"
    } else if(name_spp %in% RS$s_13_TenSugi){
      name_RS <- "高知天然スギ"
    } else if(name_spp %in% RS$s_13_Hinoki){
      name_RS <- "高知ヒノキ"
    } else if(name_spp %in% RS$s_13_TenHinoki){
      name_RS <- "高知天然ヒノキ"
    } else if(name_spp %in% RS$s_13_Akamatsu){
      name_RS <- "高知アカマツ"
    } else if(name_spp %in% RS$s_13_Momi){
      name_RS <- "高知モミ"
    } else if(name_spp %in% RS$s_13_Tsuga){
      name_RS <- "高知ツガ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "名古屋カラマツ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "高知モミ"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "高知広葉樹"
    } else {
      stop("No such species in this region.")
    }
    ## 14. Kumamoto
  } else if (name_region %in% RS$r_14_Kumamoto){
    if(name_spp %in% RS$s_14_Sugi & name_region == "飫肥"){
      name_RS <- "飫肥スギ"
    } else if(name_spp %in% RS$s_14_Sugi){
      name_RS <- "熊本スギ"
    } else if(name_spp %in% RS$s_14_TenSugi & name_region == "下屋久"){
      name_RS <- "下屋久天然スギ"
    } else if(name_spp %in% RS$s_14_TenSugi){
      name_RS <- "上屋久天然スギ"
    } else if(name_spp %in% RS$s_14_Hinoki){
      name_RS <- "熊本ヒノキ"
    } else if(name_spp %in% RS$s_14_Akamatsu & name_region == "霧島"){
      name_RS <- "霧島アカマツ"
    } else if(name_spp %in% RS$s_14_Akamatsu){
      name_RS <- "熊本アカマツ"
    } else if(name_spp %in% RS$s_14_TenAkamatsu){
      name_RS <- "霧島天然アカマツ"
    } else if(name_spp %in% RS$s_14_Momi){
      name_RS <- "熊本モミ"
    } else if(name_spp %in% RS$s_14_Tsuga){
      name_RS <- "熊本ツガ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "名古屋カラマツ"
    } else if(name_spp %in% c("イチョウ", "ナギ")){ # to include as Broadleaf
      name_RS <- "熊本広葉樹1類"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "熊本モミ"
    } else if(any(str_detect(name_spp, RS$s_14_Broad_1))){
      name_RS <- "熊本広葉樹1類"
    } else if(any(str_detect(name_spp, RS$s_14_Broad_2))){
      name_RS <- "熊本広葉樹2類"
    } else {
      stop("No such species in this region.")
    }
  } else {
    if(is.null(name_invalid)){
      stop("No such region.")
    } else {
      name_RS <- name_invalid
    }

  }

  return(name_RS)
}


#' Get Region name
#'
#' Wrapper function to get Region name that is used in [stemVolume()] (e.g., "東京スギ").
#' This is faster when large amount of data are passed and those are duplicated.
#'
#' Species & Region names are defined following Hosoda et al. (2010) Jpn. J. For. Plann. 44:23-39.
#'  (https://doi.org/10.20659/jjfp.44.2_23) and
#'  "幹材積計算プログラム" developed by FFPRI.
#' (https://www.ffpri.affrc.go.jp/database/stemvolume/index.html)
#'
#' @param Region  Character vector. Region/Prefecture name that are listed in `RS`.
#'   Please see the original article by Hosoda et al. (2010) and description in the above program.
#' @param Spp  Character vector. Species name in Japanese.
#'   If it does not match any species in the lists, then `広葉樹` is returned.
#'   Please see the original article by Hosoda et al. (2010) and description in the above program.
#' @param RS  A list of parameters, which is generated by [getRegionName()].
#'   list names should be strictly the same as one generated by [getRegionName()] to derive the correct parameters.
#' @return Name for identifying and calculating stem volume in [stemVolume()]
#'
#' @importFrom purrr pmap
#' @importFrom dplyr %>% filter mutate group_by summarize
#' @export
#'
volumeName <- function(Region, Spp, RS = NULL){
  # for single value ---
  if(length(Region) == 1 & length(Spp) == 1){
    Name <- volumeNameSingle(Region, Spp, RS)

  # for multiple values ---
  } else {
    ## get parameters for Region & Name ---
    if(is.null(RS)) RS <- getRegionName()
    ## unique combinations of Region & Spp ---
    dt_comb <- data.frame(Region = Region, Spp = Spp) |>
      group_by(Region, Spp) |>
      summarize(.groups = "drop")
    ## get Name for each unique comb. ---
    Name_comb <- pmap(dt_comb, volumeNameSingle, name_invalid = NULL)
    ## allocate to each Region & Spp ---
    Name <- factor(paste0(Region, Spp),
                   levels = paste0(dt_comb$Region, dt_comb$Spp),
                   labels = unlist(Name_comb)) |>
      as.character()
  }
  return(Name)
}




#' Get region & name as list
#'
#' @return A list of parameters for [volumeName()]
#'
getRegionName <- function(){
  # list of species
  s_x_Todomatsu <- c("トドマツ", "アカトドマツ", "アオトドマツ") # added
  s_x_Ezomatsu  <- c("エゾマツ", "アカエゾマツ", "アカエゾ", "クロエゾマツ", "クロエゾ") # added
  s_x_Karamatsu <- c("カラマツ")
  s_x_Hiba      <- c("ヒバ", "ヒノキアスナロ", "アスナロ") # added
  s_x_Buna      <- c("ブナ")
  s_x_Sugi      <- c("スギ")
  s_x_TenSugi   <- c("天然スギ", "天スギ")
  s_x_Hinoki    <- c("ヒノキ")
  s_x_TenHinoki <- c("天然ヒノキ", "天ヒノキ")
  s_x_Akamatsu  <- c("アカマツ", "クロマツ")
  s_x_ConiferOther   <- c("針葉樹", "その他針葉樹", "針", "針葉", "N",
                          "イチョウ", "(.+|)マツ", "(.+|)ゴヨウ", "(.+|)モミ", "(.+|)シラビソ", "(.+|)シラベ",
                          "(.+|)トウヒ", "(.+|)マツハダ", "ツガ", "コメツガ", "トガサワラ", "(.+|)スギ",
                          "コウヤマキ", "(.+|)ビャクシン", "ネズミサシ", "(.+|)ネズ", "ヒノキ", "サワラ",
                          "ネズコ", "クロベ", "アスナロ", "ヒノキアスナロ", "ヒバ", "アテ", "イヌマキ",
                          "ナギ", "イヌガヤ", "ハイイヌガヤ", "イチイ", "カヤ")
  s_x_BroadleafOther <- c("広葉樹", "その他広葉樹", "広", "広葉", "|.+") # match all

  list_RS <- list(
    # list of each region
    r_1_Asahikawa = c("旭川", "宗谷支庁", "留萌支庁", "上川支庁", "雨竜地方", "宗谷", "留萌", "上川北部", "上川南部"),
    r_2_Kitami = c("北見", "網走支庁", "網走西部", "網走東部"),
    r_3_Obihiro = c("帯広", "根室支庁", "釧路支庁", "十勝支庁", "十勝", "釧路根室"),
    r_4_Sapporo = c("札幌", "日高支庁", "空知支庁", "石狩支庁", "幌別地方", "白老地方", "勇払地方", "積丹地方",
                     "古平地方", "余市地方", "石狩空知", "日高"),
    r_5_Hakodate= c("函館", "胆振支庁", "後志支庁", "桧山支庁", "渡島支庁", "渡島檜山", "後志胆振", "胆振東部"),
    r_6_Aomori  = c("青森", "岩手", "宮城"),
    r_7_Akita = c("秋田", "山形"),
    r_8_Maebashi = c("前橋", "会津新潟", "会津地方", "関東中部", "福島", "栃木", "群馬", "新潟"),
    r_9_Tokyo = c("茨城", "埼玉", "千葉", "東京", "神奈川", "山梨", "静岡"),
    r_10_Nagano = c("長野"),
    r_11_Nagoya = c("名古屋", "富山", "岐阜", "愛知"),
    r_12_Osaka  = c("山陰", "石川", "福井", "三重", "滋賀", "京都", "大阪", "兵庫", "奈良", "和歌山", "鳥取",
                     "島根", "岡山", "広島", "山口",
                     "中丹地方", "丹後地方", "但馬地方"),
    r_13_Kochi = c("徳島", "香川", "愛媛", "高知"),
    r_14_Kumamoto = c("福岡", "佐賀", "長崎", "熊本", "大分", "宮崎", "鹿児島", "沖縄",
                       "飫肥", "下屋久", "上屋久", "霧島", "飫肥地方", "下屋久地方", "上屋久地方", "霧島地方"),

    # list of species
    s_x_Todomatsu = s_x_Todomatsu,
    s_x_Ezomatsu  = s_x_Ezomatsu,
    s_x_Karamatsu = s_x_Karamatsu,
    s_x_Hiba      = s_x_Hiba,
    s_x_Buna      = s_x_Buna,
    s_x_Sugi      = s_x_Sugi,
    s_x_TenSugi   = s_x_TenSugi,
    s_x_Hinoki    = s_x_Hinoki,
    s_x_TenHinoki = s_x_TenHinoki,
    s_x_Akamatsu  = s_x_Akamatsu,
    s_x_ConiferOther    = s_x_ConiferOther,
    s_x_BroadleafOther  = s_x_BroadleafOther,
    s_1_Todomatsu = s_x_Todomatsu,
    s_1_Ezomatsu  = s_x_Ezomatsu,
    s_2_Todomatsu = s_x_Todomatsu,
    s_2_Ezomatsu  = s_x_Ezomatsu,
    s_3_Todomatsu = s_x_Todomatsu,
    s_3_Ezomatsu  = s_x_Ezomatsu,
    s_4_Todomatsu = s_x_Todomatsu,
    s_4_Ezomatsu  = s_x_Ezomatsu,
    s_5_Todomatsu = s_x_Todomatsu,
    s_5_Ezomatsu  = s_x_Ezomatsu,
    s_5_Hiba      = s_x_Hiba,
    s_5_Buna      = s_x_Buna,
    s_6_Sugi      = s_x_Sugi,
    s_6_TenSugi   = s_x_TenSugi,
    s_6_Akamatsu  = s_x_Akamatsu,
    s_6_Hiba      = s_x_Hiba,
    s_7_Sugi      = s_x_Sugi,
    s_7_TenSugi   = s_x_TenSugi,
    s_7_Karamatsu = s_x_Karamatsu,
    s_7_Akamatsu  = s_x_Akamatsu,
    s_8_Sugi      = s_x_Sugi,
    s_8_Hinoki    = s_x_Hinoki,
    s_8_Karamatsu = s_x_Karamatsu,
    s_8_Akamatsu  = s_x_Akamatsu,
    s_k_Momi      = c("モミ", "ウラジロモミ", "シラベ", "トウヒ"),
    s_k_Tsuga     = c("ツガ", "コメツガ", "ネズコ"),
    s_9_Sugi      = s_x_Sugi,
    s_9_Hinoki    = s_x_Hinoki,
    s_9_Akamatsu  = c(s_x_Akamatsu, "リュウキュウマツ"),
    s_10_Sugi      = s_x_Sugi,
    s_10_Hinoki    = s_x_Hinoki,
    s_10_TenHinoki = s_x_TenHinoki,
    s_10_Sawara    = c("サワラ", "ヒバ", "ネズコ"),
    s_10_Akamatsu  = s_x_Akamatsu,
    s_10_Karamatsu = s_x_Karamatsu,
    s_10_Buna      = c("ブナ", "ナラ", "カンバ"),
    s_10_Sawagurumi= c("サワグルミ", "ホオノキ", "カツラ", "シオジ"),
    s_11_Sugi      = s_x_Sugi,
    s_11_Hinoki    = s_x_Hinoki,
    s_11_TenHinoki = s_x_TenHinoki,
    s_11_Sawara    = c("サワラ", "ヒバ"),
    s_11_Akamatsu  = s_x_Akamatsu,
    s_11_Karamatsu = s_x_Karamatsu,
    s_12_Sugi      = s_x_Sugi,
    s_12_TenSugi   = s_x_TenSugi,
    s_12_Hinoki    = s_x_Hinoki,
    s_12_Akamatsu  = s_x_Akamatsu,
    s_12_Kouyamaki = c("コウヤマキ", "サワラ"),
    s_12_Momi      = c("モミ"),
    s_12_Broad_1   = c("サワグルミ", "オニグルミ", "ホオノキ", "カツラ", "クヌギ", "アベマキ", "センノキ", "シオジ",
                        "ミズキ", "キハダ", "トネリコ", "アサダ", "ヤチダモ", "ニレ", "キリ", "ドロノキ"),
    s_12_Broad_2   = s_x_BroadleafOther,
    s_13_Sugi      = s_x_Sugi,
    s_13_TenSugi   = s_x_TenSugi,
    s_13_Hinoki    = s_x_Hinoki,
    s_13_TenHinoki = s_x_TenHinoki,
    s_13_Akamatsu  = s_x_Akamatsu,
    s_13_Momi      = c("モミ"),
    s_13_Tsuga     = c("ツガ"),
    s_14_Sugi      = s_x_Sugi,
    s_14_TenSugi   = s_x_TenSugi,
    s_14_Hinoki    = s_x_Hinoki,
    s_14_Akamatsu  = s_x_Akamatsu,
    s_14_TenAkamatsu=c("天然アカマツ"),
    s_14_Momi      = c("モミ"),
    s_14_Tsuga     = c("ツガ"),
    s_14_Broad_1   = c(c("サワグルミ", "シデ", "クリ", "クヌギ", "シイ", "ケヤキ", "カツラ", "ホオノキ", "エンジュ",
                          "センダン", "アブラギリ", "イイギリ", "ハリギリ", "ミヤコダラ", "ミズキ", "ヤマガキ", "トネリコ",
                          "シオジ", "チシャノキ", "ナギ", "イチョウ"),
                        c(".+シイ", ".+ジイ")),
    s_14_Broad_2   = s_x_BroadleafOther
  )
  return(list_RS)
}





#' Generate a list of dataframe for parameters
#'
#' Generate dataframe of parameters for equations 1) stem volume, 2) BH form factor, and
#' 3) height form based volume calculation.
#'
#' @return A list of dataframe for [stemVolume()]
#'
#' @importFrom dplyr %>% filter mutate
#' @export
#'
getStemCoefficients <- function(){
  # coefficients for stem volume calculation ---
  dt_stem <- data.frame(
    EW = c('東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東', '東',
           '東', '東', '東', '東', '東', '東', '東', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西', '西',
           '西', '西', '西', '西', '西', '西', '西', '西', '西', '西'),
    region = c('旭川', '旭川', '旭川', '旭川', '旭川', '旭川', '北見', '北見', '北見', '北見', '北見', '北見', '北見',
               '北見', '帯広', '帯広', '帯広', '帯広', '帯広', '帯広', '帯広', '帯広', '帯広', '札幌', '札幌', '札幌',
               '札幌', '札幌', '札幌', '函館', '函館', '函館', '函館', '函館', '函館', '函館', '函館', '函館', '函館',
               '函館', '函館', '函館', '函館', '北海道', '北海道', '北海道', '北海道', '青森', '青森',
               '青森', '青森', '青森', '青森', '青森', '青森', '青森', '青森', '青森', '秋田', '秋田', '秋田', '秋田',
               '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '秋田', '関東中部',
               '関東中部', '関東中部', '関東中部', '関東中部', '関東中部', '関東中部', '関東中部', '関東中部', '関東中部',
               '関東中部', '前橋', '前橋', '前橋', '前橋', '前橋', '会津新潟', '会津新潟', '会津新潟', '前橋', '前橋',
               '前橋', '前橋', '前橋', '前橋', '前橋', '前橋', '会津新潟', '会津新潟', '前橋', '前橋', '前橋', '前橋',
               '前橋', '前橋', '前橋', '前橋', '東京', '東京', '東京', '東京', '東京', '東京', '東京', '東京', '東京',
               '東京', '東京', '東京', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野',
               '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野',
               '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '長野', '名古屋', '名古屋', '名古屋',
               '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋',
               '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋', '名古屋',
               '名古屋', '名古屋', '大阪', '山陰', '大阪', '大阪', '大阪', '大阪', '山陰', '山陰', '山陰', '大阪', '大阪',
               '大阪', '大阪', '大阪', '大阪', '大阪', '大阪', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知',
               '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知',
               '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '高知', '熊本',
               '熊本', '飫肥', '飫肥', '飫肥', '飫肥', '飫肥', '下屋久', '下屋久', '下屋久', '下屋久', '下屋久', '熊本',
               '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '霧島', '霧島',
               '霧島', '霧島', '霧島', '霧島', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本', '熊本'),
    species = c('トドマツ', 'トドマツ', 'トドマツ', 'トドマツ', 'トドマツ', 'エゾマツ', 'トドマツ', 'トドマツ', 'トドマツ',
                'エゾマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'トドマツ', 'トドマツ', 'トドマツ', 'トドマツ',
                'トドマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'トドマツ', 'トドマツ', 'トドマツ', 'トドマツ',
                'エゾマツ', 'エゾマツ', 'トドマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'エゾマツ', 'ヒバ',
                'ヒバ', 'ヒバ', 'ヒバ', 'ブナ', 'ブナ', 'ブナ', 'ブナ', '広葉樹', '広葉樹', '広葉樹', '広葉樹', 'スギ',
                'スギ', 'スギ', 'スギ', 'アカマツ', 'アカマツ', 'アカマツ', '広葉樹', '広葉樹', '広葉樹', '広葉樹', 'スギ',
                'スギ', 'スギ', 'スギ', '天然スギ', '天然スギ', '天然スギ', '天然スギ', 'カラマツ', 'カラマツ', 'カラマツ',
                'カラマツ', 'カラマツ', '広葉樹', '広葉樹', '広葉樹', 'モミ', 'モミ', 'モミ', 'モミ', 'モミ', 'ツガ', 'ツガ',
                'ツガ', 'ツガ', 'ツガ', 'ツガ', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ', 'ヒノキ',
                'ヒノキ', 'ヒノキ', 'ヒノキ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ',
                'カラマツ', 'カラマツ', 'カラマツ', 'カラマツ', 'カラマツ', '広葉樹', '広葉樹', '広葉樹', 'スギ', 'スギ',
                'スギ', 'スギ', 'ヒノキ', 'ヒノキ', 'ヒノキ', 'アカマツ', 'アカマツ', 'アカマツ', '広葉樹', '広葉樹', 'スギ',
                'スギ', 'スギ', 'スギ', 'ヒノキ', 'ヒノキ', 'ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ',
                '天然ヒノキ', '天然ヒノキ', 'サワラ', 'サワラ', 'サワラ', 'サワラ', 'アカマツ', 'アカマツ', 'アカマツ',
                'アカマツ', 'カラマツ', 'カラマツ', 'カラマツ', 'カラマツ', 'ブナ', 'ブナ', 'ブナ', 'ブナ', 'サワグルミ',
                'サワグルミ', 'サワグルミ', 'サワグルミ', 'スギ', 'スギ', 'スギ', 'スギ', 'ヒノキ', 'ヒノキ', 'ヒノキ',
                'ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', 'サワラ', 'サワラ', 'アカマツ', 'アカマツ',
                'アカマツ', 'カラマツ', 'カラマツ', 'カラマツ', 'カラマツ', 'カラマツ', '広葉樹', '広葉樹', '広葉樹',
                '広葉樹', '広葉樹', 'スギ', '天然スギ', 'ヒノキ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ',
                'アカマツ', 'コウヤマキ', 'モミ', '広葉樹1型', '広葉樹1型', '広葉樹1型', '広葉樹2型', '広葉樹2型',
                '広葉樹2型', 'スギ', 'スギ', 'スギ', 'スギ', '天然スギ', '天然スギ', '天然スギ', 'ヒノキ', 'ヒノキ', 'ヒノキ',
                'ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', '天然ヒノキ', 'アカマツ',
                'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'モミ', 'モミ', 'モミ', 'ツガ', 'ツガ', '広葉樹',
                '広葉樹', '広葉樹', '広葉樹', '広葉樹', '広葉樹', '広葉樹', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ', 'スギ',
                'スギ', '天然スギ', '天然スギ', '天然スギ', '天然スギ', '天然スギ', 'ヒノキ', 'ヒノキ', 'ヒノキ', 'モミ',
                'モミ', 'モミ', 'ツガ', 'ツガ', 'ツガ', 'ツガ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ', 'アカマツ',
                'アカマツ', '天然アカマツ', '天然アカマツ', '広葉樹1類', '広葉樹1類', '広葉樹2類', '広葉樹2類',
                '広葉樹2類', '広葉樹2類', '広葉樹2類', '広葉樹2類'),
    var_left = c('log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'v', 'v', 'v', 'v', 'v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'v', 'v', 'v', 'v', 'log v', 'log v', 'log v', 'v',
                 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'v', 'v', 'v', 'v', 'v', 'v', 'v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v',
                 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v', 'log v'),
    D_lower = c(0, 21, 31, 41, 51, 0, 0, 11, 31, 0, 11, 21, 31, 71, 0, 11, 21, 41, 51, 0, 31, 41, 51, 0, 11, 21, 51, 0, 31,
                0, 0, 51, 61, 71, 81, 0, 11, 21, 41, 0, 11, 31, 51, 0, 11, 21, 31, 0, 11, 21, 41, 0, 23, 47, 0, 11, 51, 71,
                0, 11, 21, 41, 0, 11, 21, 71, 0, 11, 21, 31, 41, 0, 11, 45, 0, 11, 21, 31, 51, 0, 11, 21, 31, 41, 51, 0, 11,
                21, 31, 41, 0, 11, 31, 0, 11, 21, 31, 0, 11, 21, 41, 0, 41, 0, 11, 21, 31, 41, 0, 11, 41, 0, 11, 31, 41, 0,
                11, 21, 0, 31, 41, 0, 61, 0, 11, 21, 31, 0, 11, 21, 0, 11, 21, 31, 41, 51, 0, 11, 21, 31, 0, 11, 21, 41, 0,
                11, 21, 31, 0, 21, 31, 41, 0, 11, 21, 31, 0, 11, 31, 41, 0, 11, 21, 41, 0, 21, 31, 41, 0, 51, 0, 11, 21, 0,
                11, 21, 31, 41, 0, 11, 21, 31, 61, 0, 0, 0, 0, 11, 31, 0, 11, 41, 0, 0, 0, 11, 21, 0, 11, 21, 0, 11, 21, 31,
                0, 41, 101, 0, 11, 21, 31, 0, 21, 31, 41, 51, 61, 0, 11, 21, 31, 41, 51, 0, 51, 81, 0, 81, 0, 11, 21, 31, 41,
                51, 61, 0, 31, 0, 11, 21, 41, 51, 0, 31, 71, 81, 91, 0, 11, 21, 0, 21, 91, 0, 11, 21, 71, 0, 21, 0, 11, 21,
                31, 0, 41, 0, 51, 0, 11, 21, 31, 41, 51),
    D_upper = c(21, 31, 41, 51, 999, 999, 11, 31, 999, 11, 21, 31, 71, 999, 11, 21, 41, 51, 999, 31, 41, 51, 999, 11, 21, 51,
                999, 31, 999, 999, 51, 61, 71, 81, 999, 11, 21, 41, 999, 11, 31, 51, 999, 11, 21, 31, 999, 11, 21, 41, 999,
                23, 47, 999, 11, 51, 71, 999, 11, 21, 41, 999, 11, 21, 71, 999, 11, 21, 31, 41, 999, 11, 45, 999, 11, 21, 31,
                51, 999, 11, 21, 31, 41, 51, 999, 11, 21, 31, 41, 999, 11, 31, 999, 11, 21, 31, 999, 11, 21, 41, 999, 41, 999,
                11, 21, 31, 41, 999, 11, 41, 999, 11, 31, 41, 999, 11, 21, 999, 31, 41, 999, 61, 999, 11, 21, 31, 999, 11, 21,
                999, 11, 21, 31, 41, 51, 999, 11, 21, 31, 999, 11, 21, 41, 999, 11, 21, 31, 999, 21, 31, 41, 999, 11, 21, 31,
                999, 11, 31, 41, 999, 11, 21, 41, 999, 21, 31, 41, 999, 51, 999, 11, 21, 999, 11, 21, 31, 41, 999, 11, 21, 31,
                61, 999, 999, 999, 999, 11, 31, 999, 11, 41, 999, 999, 999, 11, 21, 999, 11, 21, 999, 11, 21, 31, 999, 41, 101,
                999, 11, 21, 31, 999, 21, 31, 41, 51, 61, 999, 11, 21, 31, 41, 51, 999, 51, 81, 999, 81, 999, 11, 21, 31, 41,
                51, 61, 999, 31, 999, 11, 21, 41, 51, 999, 31, 71, 81, 91, 999, 11, 21, 999, 21, 91, 999, 11, 21, 71, 999, 21,
                999, 11, 21, 31, 999, 41, 999, 51, 999, 11, 21, 31, 41, 51, 999),
    coef_a = c(-4.13703479312368, -4.44880936868777, -4.1996051877536, -4.34006490330827, -4.00382390330827, -4.1502711877536,
               -4.10644380865396, -4.1619207695884, -4.12663926964221, -4.05552967471382, -4.0860581877536, -4.20052936868777,
               -4.22264167471382, -3.61188382617647, -3.95894879312368, -4.24844452944638, -4.00200965914322, -4.00219052731017,
               -3.87760811910794, -4.05920537189881, -4.07402684117393, -4.49781401928009, -4.10699215377484, -4.03879629972297,
               -4.10528447651336, -4.13156426964221, -3.67793360159354, -4.21723326964221, -4.0185289850586, -4.09519329595912,
               -4.07232037360688, -3.48443653730748, -3.98045005891573, -3.307466287191, -3.90972037360688, -4.0084211,
               -4.1072282, -4.2521633, -3.979531, -4.13888374, -4.3324403, -4.2630919, -4.2585898, -4.068644, -4.335395,
               -4.441199, -4.332596, -4.175185, -4.270204, -4.205519, -3.854066, -4.233, -4.23317, -4.86886, -4.1916, -4.36445,
               -4.11234, -3.93414, -4.117135, -4.221487, -4.326722, -4.072908, -4.0807061, -4.2395022, -4.1558594, -3.8600114,
               NA, NA, NA, NA, NA, -4.1627841, -4.3497673, -4.1142834, -4.094453, -4.151762, -4.200651, -4.192596, -4.250481,
               -4.25082, -4.358637, -4.338446, -4.112924, -4.073435, -4.090918, -4.1231, -4.26496, -4.28486, -4.17044,
               -4.11774, -4.0989, -4.27591, -4.09301, -4.173533, -4.293729, -4.271259, -4.404407, -4.249503, -4.155639,
               -4.194535, -4.42347, -4.144592, -4.356181, -4.155099, -4.369281, -4.348104, -3.976731, -4.431495, -4.20067295,
               -4.32216295, -4.15096808,  -4.172632, -4.219069, -4.211821, -3.921218, -4.230151, -4.354186, -4.229853,
               -4.249808, -4.060353, -4.347438, -4.344385, -4.174703, -4.210459, -4.227928, -4.321558, -4.185437,
               -4.12648327360688, -4.3680436850586, -4.19750910159354, -4.2263829, -4.3714161, -4.4342693, -4.3824548,
               -3.9865606, -4.106319, -4.2091652, -4.1875621, -4.7386515, -4.2061745, -4.2361, -4.28665, -4.34632, -3.93698,
               -4.22457229972296, -4.4138789850586, -4.32742915906381, -4.20859568413365, -4.26067399774047,
               -4.41706146789049, -4.60811372701199, -4.22251315479781, -4.27548592701199, -4.34448094735743,
               -4.50207783292256, -4.44378081928009, -4.229266, -4.265222, -4.068185, -3.923549, -4.341504, -4.416053,
               -4.270848, -3.75679, -4.2968, -4.0799, -4.39928, -4.13316, -4.185938, -4.028179, -4.333173, -4.258013,
               -4.139304, -4.217893, -4.391799, -4.36818, -4.417533, -4.187525, -4.166839, -4.366075, -4.532942, -4.186608,
               -3.785556, -4.19207, -4.14295, -4.31101, -4.14929, -4.10009, -4.342135, -4.134795, -4.097423,
               -4.202034, -4.16774, -4.1576, -4.122284, -4.28339, -4.418992, -4.070481, -4.232323, -4.272709, 6.901e-05,
               0.00177071, 0.00442833, 0.02464082, -4.07379531774047, -4.00703928185152, -3.01596792824349, 0.00016148,
               0.00044496, 0.00526351, 0.01054664, 0.00201923, 0.00808959, -0.01648444, 0.04483989, 0.05748548, 0.02273942,
               -0.00028, 0.000558, 0.001466, 0.004393, 0.004158, -0.005548, -4.22931629985458, -4.21578030502934,
               -4.04234599383693, -4.17986612070474, -3.94444833261126, 9.1e-05, 0.00185, 0.002574, 0.02873, -0.025002,
               -0.006166, 0.045811, -4.203818, -3.9245239, -4.1569, -4.175079, -4.005031, -4.253183, -4.6572346, -4.114226,
               -4.2231395, -4.0174538, -3.9272498, -3.092297, -4.12789, -4.317069, -4.2014653, -4.171861, -4.2536036,
               -3.468589, -4.14579, -4.425594, -4.138925, -3.8941264, -4.1209, -4.1120801, -4.11513, -4.245534, -4.381961,
               -4.1340831, -4.894381, -4.2002787, -4.1992744, -3.7806128, -4.06484, -4.224395, -4.080059, -4.178527,
               -4.046078, -3.9095253),
    var_a = c('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', 'h', 'h', 'h', 'h', '', '', '', 'h', 'h', 'h', 'h', 'h', 'h',
              'h', 'h', 'h', 'h', 'h', 'h', 'h', 'h', 'h', 'h', '', '', '', '', '', 'h', 'h', 'h', 'h', 'h', 'h', 'h',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', ''),
    coef_b = c(1.93071, 2.048569, 1.808385, 1.908367, 1.641559, 1.796738, 1.90563, 1.865311, 1.778222, 1.954285, 1.780146,
               1.864287, 1.827392, 1.409974, 1.63873, 1.970586, 1.947286, 1.705271, 1.711802, 1.820048, 2.076056, 2.172937,
               1.776975, 1.81669, 1.776303, 1.797374, 1.594803, 1.873659, 1.79369, 1.681121, 1.82408, 1.568947, 1.783992,
               1.495464, 1.723599, 1.8054881, 1.8483473, 1.7922954, 1.6779623, 1.7674314, 1.9039982, 1.8162864, 1.8718783,
               1.756152, 1.903051, 1.853014, 1.848675, 1.930072, 1.933985, 1.778137, 1.616131, 1.93399, 1.75974, 1.802,
               1.84673, 1.87135, 1.6384, 1.75774, 1.769161, 1.810503, 1.726305, 1.617248, 1.8757208, 1.8706824, 1.8234115,
               1.7086843, NA, NA, NA, NA, NA, 0.9096719, 0.9715385, 0.9200743, 1.953184, 1.922623, 1.931306, 1.858442,
               1.816918, 2.04338, 2.087948, 2.025199, 1.833615, 1.825181, 1.701661, 1.7454, 1.83346, 1.74357, 1.76381,
               1.79204, 1.8876, 1.82147, 1.63363, 1.810672, 1.905709, 1.771888, 1.671519, 1.946501, 1.847898, 1.804172,
               2.006485, 1.867682, 1.91579, 1.80389, 1.979213, 1.904401, 1.640825, 1.90887, 1.94019664, 1.93813902,
               1.82464098, 1.753904, 1.849344, 1.944187, 1.600066, 1.906941, 1.867405, 1.710344, 1.863288, 1.895653,
               2.002385, 1.834463, 1.765531, 1.83883, 1.84978, 1.842796, 1.703437, 1.8168788, 1.8121862, 1.6603092, 1.9709522,
               1.968849, 1.9418846, 1.9125607, 1.7022045, 1.7516722, 1.9183687, 1.7120116, 2.1613151, 1.7750576, 1.77,
               1.78364, 1.88856, 1.683851, 1.87385, 1.96416, 1.84173, 1.74034, 1.8447874, 1.8292605, 1.8435209, 1.776932,
               1.7783579, 1.8069719, 1.891904, 1.8204232, 1.967735, 1.864665, 1.687367, 1.728859, 2.009104, 1.909076, 1.7564,
               1.370193, 1.88715, 1.7777, 1.95151, 1.74574, 1.864437, 1.78424, 1.969025, 1.880652, 1.793262, 1.839235,
               1.995252, 1.873209, 1.688448, 1.646844, 1.852021, 1.857805, 2.048645, 1.846988, 1.585865, 1.82696, 1.95677,
               1.83546, 1.980643, 1.814939, 1.920082, 1.835981, 1.871126, 1.838879, 1.75216, 1.86611, 1.744307, 1.749564,
               1.681823, 1.856641, 1.864235, 1.752091, 0.58810351, 1.04476089, 1.42509179, 4.65164113, 1.78647567, 1.7644211,
               1.14155393, 0.28229011, 0.1102957, 0.70722739, 0.5264662, 0.35600997, 1.73146183, -2.2243528, 3.65830858,
               3.61919679, 0.30490816, NA, NA, NA, NA, NA, NA, 1.7641027, 1.6984878, 1.6417125, 1.84841224, 1.68040853,
               0.208874, 1.190725, 0.216124, 3.364846, -3.488326, 0.899678, 1.599526, 1.819629, 1.6644424, 1.95461, 1.908431,
               1.730365, 1.7833516, 2.0270131, 1.932645, 1.9620566, 1.7995275, 1.7838516, 1.4551198, 1.93699, 1.921617,
               1.786204, 2.000873, 1.8556134, 1.4779742, 2.06848, 2.119952, 1.9579086, 1.7406157, 1.881206, 1.8825315,
               1.90971, 1.820157, 1.781114, 1.6292859, 1.986659, 1.8696292, 1.858214, 1.7774304, 1.86248, 1.883197,
               1.769327, 1.85037, 1.7631913, 1.6181613),
    var_b = c('log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', '', '', '', '', '', 'd2h', 'd2h', 'd2h', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'd2', 'd2', 'd2', 'd2', 'log d', 'log d', 'log d', 'd2', 'd2',
              'd2', 'd2', 'd2', 'd2', 'd2', 'd2', 'd2', 'd2', '', '', '', '', '', '', 'log d', 'log d', 'log d', 'log d',
              'log d', 'd2', 'd2', 'd2', 'd2', 'd2', 'd2', 'd2', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d',
              'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d', 'log d'),
    coef_c = c(0.884533, 0.992988, 1.06837, 1.045786, 1.126041, 1.032912, 0.870334, 0.976962, 1.051595, 0.732053, 0.98618,
               0.997385, 1.0565, 1.166567, 0.979937, 0.928967, 0.741914, 1.025535, 0.950097, 0.904747, 0.635918, 0.82625,
               1.024414, 0.934425, 1.038109, 1.03448, 0.964676, 1.019256, 0.956599, 1.131348, 0.934568, 0.842787, 0.931472,
               0.845897, 0.967339, 0.8620075, 0.9269954, 1.1063961, 1.0363396, 0.963137, 1.005558, 1.0610156, 0.986726,
               0.90621, 1.02541, 1.166956, 1.088954, 0.890582, 0.974995, 1.086571, 1.026974, 0.92679, 1.12264, 1.52925,
               0.92548, 1.06862, 1.17842, 0.903, 0.97415, 1.044206, 1.227196, 1.170206, 0.8353226, 1.0122438, 1.0059217,
               0.951969, 0.37872483, 0.37777926, 0.35967623, 0.36363858, 0.33591133, NA, NA, NA, 0.784224, 0.883683,
               0.91936, 1.00008, 1.095799, 0.836823, 0.873856, 0.915023, 0.957762, 0.940891, 1.110799, 1.014, 1.06569,
               1.17719, 1.06412, 0.99303, 0.8223, 1.0877, 1.16324, 0.982833, 1.011385, 1.138415, 1.363617, 0.942682, 0.951955,
               1.034248, 0.967757, 0.930308, 1.02159, 0.962587, 0.998347, 1.062478, 1.080387, 1.088002, 0.84689666,
               0.96697002, 0.97625989, 1.040853, 1.008086, 0.894801, 1.075361, 0.942734, 1.108487, 1.175119, 1.004738,
               0.811988, 0.888616, 1.098828, 1.073801, 1.000611, 1.010122, 1.092369, 1.151114, 0.9167444, 1.1783766,
               1.2042435, 0.8801095, 1.016603, 1.0942478, 1.081489, 1.0390683, 1.0694202, 0.8966176, 1.1141249, 1.0727478,
               1.0946404, 1.0987, 1.13056, 1.067932, 0.985842, 0.94852, 1.04523, 1.1108, 1.13316, 1.0088782, 1.1567557,
               1.2867129, 1.0830437, 1.0725259, 1.1166818, 1.1516396, 1.1883403, 0.874649, 1.023757, 1.079349, 0.927572,
               0.929831, 1.097628, 1.166107, 1.239011, 1.0419, 0.99996, 1.04225, 1.09378, 0.977728, 0.96031, 0.989459,
               0.990775, 0.995431, 0.961745, 0.986175, 1.094553, 1.334207, 1.226351, 0.896175, 1.084483, 1.013891, 0.972756,
               1.025991, 0.99227, 0.82008, 1.10655, 0.811703, 0.91831, 0.989263, 0.939368, 0.863667, 0.98589, 1.0711,
               0.97264, 0.962038, 1.129285, 1.309814, 0.819044, 0.973986, 1.131128, 0.38337273, 0.29644026, 0.25603657,
               0.11705915, 0.95616277, 0.9422713, 1.12285567, 0.37512815, 0.37336909, 0.2961912, 0.27740293, 0.33508684,
               0.23774503, 0.47430445, 0.13666666, 0.14882345, 0.26804944, 0.435094, 0.344742, 0.321958, 0.298333, 0.295753,
               0.328098, 1.1414065, 1.2200706, 1.1826469, 0.9925352, 1.04814308, 0.371003, 0.25383, 0.338284, 0.144179,
               0.43199, 0.292977, 0.220153, 1.025738, 0.9881512, 0.84007, 0.910503, 0.966908, 1.0895793, 1.0710697, 0.831036,
               0.8723995, 0.9407524, 0.8962587, 0.7782647, 0.81243, 1.016795, 1.0696647, 0.833059, 1.0611038, 1.0397065,
               0.67053, 0.896909, 0.8223236, 0.9457805, 0.885486, 0.8631743, 0.80472, 1.045666, 1.21021, 1.2063186, 1.385014,
               1.0046008, 0.9851158, 0.7698188, 0.83054, 0.967565, 0.973818, 0.971142, 0.9778835, 1.0823696),
    var_c = c('log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', '', '', '', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'd2h', 'd2h', 'd2h', 'd2h', 'log h', 'log h', 'log h', 'd2h',
              'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'd2h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h',
              'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h', 'log h'),
    coef_d = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.00323, 0.00577, 0.0212, 0.00048, 0.12648, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, -0.00018, -0.01266, -0.04404, -0.51335,
               NA, NA, NA, 0.00034, 0.00081, -0.04223, -0.08798, -0.01169, -0.09417, 0.31689, -0.79568, -1.13871, -0.07981,
               0.00257, 0.00594, 0.01261, 0.00513, 0.01454, 0.0225, NA, NA, NA, NA, NA, -3e-05, -0.01444, -0.03328, -0.48339,
               0.81153, 0.0731, -0.91675, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
               NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
    var_d = c('', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '',
              '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', ''),
    D_unit = c('cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'm', 'm', 'm', 'm', 'm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'm', 'm', 'm', 'm', 'cm', 'cm', 'cm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'm',
               'm', 'm', 'm', 'cm', 'cm', 'cm', 'cm', 'cm', 'm', 'm', 'm', 'm', 'm', 'm', 'm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm',
               'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm', 'cm')
  ) |>
    mutate(name = paste0(region, species)) %>%
    rbind(filter(., name == "下屋久天然スギ") |>
            mutate(region = "上屋久", name = "上屋久天然スギ"))


  # data for breast height form factor -----
  dt_ff <- data.frame(
    name = c('青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ',
             '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ', '青森天然スギ'),
    DBH = c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56,
            58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100, 102, 104, 106, 108,
            110, 112, 114, 116, 118, 120, 122, 124, 126, 128, 130),
    factor = c(0.578, 0.554, 0.536, 0.522, 0.51, 0.497, 0.487, 0.479, 0.472, 0.466, 0.461, 0.457,
               0.455, 0.453, 0.451, 0.45, 0.449, 0.447, 0.445, 0.441, 0.439, 0.436, 0.433, 0.43, 0.427, 0.423, 0.42, 0.417,
               0.415, 0.412, 0.41, 0.408, 0.405, 0.403, 0.4, 0.398, 0.396, 0.394, 0.391, 0.388, 0.385, 0.383, 0.381, 0.378,
               0.375, 0.373, 0.37, 0.367, 0.364, 0.362, 0.359, 0.357, 0.354, 0.351, 0.349, 0.346, 0.344, 0.341, 0.338,
               0.336, 0.333, 0.331, 0.328)
  )

  # data for height form -----
  dt_hf <- data.frame(
    name = c('青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ',
             '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ',
             '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森ヒバ',
             '青森ヒバ', '青森ヒバ', '青森ヒバ', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹',
             '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹',
             '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹',
             '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹',
             '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹', '青森針葉樹',
             '青森針葉樹', '青森針葉樹', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ', '秋田アカマツ',
             '秋田アカマツ'),
    H = c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
          36, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
          31, 32, 33, 34, 35, 36, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
          26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42),
    factor = c(3.927, 4.456, 4.941, 5.39, 5.819, 6.228, 6.617, 7, 7.365, 7.728, 8.075, 8.442, 8.778, 9.14, 9.492, 9.834,
               10.189, 10.56, 10.9, 11.284, 11.637, 12.012, 12.412, 12.78, 13.175, 13.568, 13.992, 14.382, 14.805, 15.264,
               1.9, 2.49, 3, 3.45, 3.84, 4.27, 4.72, 5.13, 5.6, 6.05, 6.48, 6.89, 7.28, 7.65, 8, 8.33, 8.64, 8.93, 9.4,
               9.87, 10.12, 10.58, 11.04, 11.25, 11.7, 12.15, 12.32, 12.76, 13.2, 13.33, 13.76, 14.19, 14.28, 14.7, 15.12,
               1.882, 2.592, 3.168, 3.6, 3.966, 4.48, 4.88, 5.301, 5.72, 6.127, 6.528, 6.916, 7.35, 7.77, 8.16, 8.568,
               8.964, 9.348, 9.72, 10.08, 10.428, 10.764, 11.064, 11.4, 11.726, 12.042, 12.348, 12.673, 12.93, 13.237,
               13.536, 13.86, 14.144, 14.49, 14.832, 15.133, 15.504, 15.834, 16.2, 16.564, 16.884)
  )

  return(
    list(stem = dt_stem,
         ff = dt_ff,
         hf = dt_hf)
  )
}


