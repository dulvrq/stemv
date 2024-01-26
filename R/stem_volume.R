
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
#' - The significant digits for the parameters of several equations are different
#'  (e.g., `函館エゾマツ`, `東京広葉樹`, `長野カラマツ`).
#' - The additional revisions on `札幌トドマツ`,  `高知天然スギ`, and `青森アカマツ` are not reflected.
#' - The range for 5 points moving average adjustment is different for `青森広葉樹` and `高知広葉樹`.
#' - The error in floating-point numbers can cause the differences in rounding calculation for `北海道針葉樹`.
#'
#' @param Name  Japanese character for identifying the calculation methods (e.g., `東京スギ`).
#'   This is expected from [volumeName()].
#' @param D  Numeric. DBH in cm.
#' @param H  Numeric. Tree height in m.
#' @param list_data  A list of data for coefficients. This is expected from [getStemCoefficients()].
#'   If NULL, internally call [getStemCoefficients()].
#' @param off_adj  Logical. If TRUE, all moving average adjustment for junctions become inactive.
#' @param ... Additional arguments.
#'
#'  These can be:
#'  * `stop_if_NA`: logical. If TRUE, cause error when NAs/NaNs/Inf/-Inf are provided for `Name`/`D`/`H`;
#'    otherwise return NAs.
#'  * `adjust_to_excel`: logical. Calculate volume in the same way as program.
#' @return Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#'
#' @importFrom dplyr %>% filter mutate group_by summarize arrange
#' @importFrom utils tail
#' @importFrom cli cli_abort cli_alert_info cli_div
#' @export
#'

stemVolumeSingle <- function(Name, D, H, list_data = NULL, off_adj = FALSE, ...){
  ## check D, H --
  cli_div(theme = list(.tmp = list(color = "yellow4", "font-style" = "italic", "font-weight" = "bold"), # set color
                       .tmp2 = list(color = "blue", "font-style" = "italic")))
  ### check length ---
  if(any(c(length(Name) != 1, length(D) != 1, length(H) != 1)))
    cli_abort(c("x" ="Only single observation is allowed for Name/D/H."))
  ### check NA/NaN/Inf ---
  if("stop_if_NA" %in% names(list(...)) && list(...)[["stop_if_NA"]] == TRUE){
    if(is.na(Name)) cli_abort(c("NA/NaN should be removed from {.tmp Name}.", "x" = "{.tmp Name} is NA/NaN."))
    if(is.na(D))    cli_abort(c("NA/NaN should be removed from {.tmp D}.", "x" = "{.tmp D} is NA/NaN."))
    if(is.na(H))    cli_abort(c("NA/NaN should be removed from {.tmp H}.", "x" = "{.tmp H} is NA/NaNs."))
    if(is.infinite(D)) cli_abort(c("Inf/-Inf should be removed from {.tmp D}.", "x" = "{.tmp D} is Inf/-Inf."))
    if(is.infinite(H)) cli_abort(c("Inf/-Inf should be removed from {.tmp H}.", "x" = "{.tmp H} is Inf/-Inf."))
  } else {
    if(is.na(Name)) cli_alert_warning("{.tmp Name} is NA/NaN.")
    if(is.na(D))    cli_alert_warning("{.tmp D} is NA/NaN.")
    if(is.na(H))    cli_alert_warning("{.tmp H} is NA/NaN.")
    if(is.infinite(D)) cli_alert_warning("{.tmp D} is Inf/-Inf.")
    if(is.infinite(H)) cli_alert_warning("T{.tmp H} is Inf/-Inf.")
  }
  ### check class --
  if(!is.numeric(D) & !is.na(D)) cli_abort(c("{.tmp D} should be numeric.", "x" = "class {.tmp D} is {.cls {class(D)}}."))
  if(!is.numeric(H) & !is.na(H)) cli_abort(c("{.tmp H} should be numeric.", "x" = "class {.tmp H} is {.cls {class(H)}}."))

  ## adjustment to D & H. if either is NA, return NA ---
  if(is.na(D) | is.na(H)) return(NA)
  if(is.infinite(D) | is.infinite(H)) return(NA)
  if(D < 0){
    cli_alert_warning("{.tmp D} is negative values (< 0), adjusting to 0.")
    D <- 0
  }
  if(H < 0){
    cli_alert_warning("{.tmp H} is negative values (< 0), adjusting to 0.")
    H <- 0
  }

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
    if(nrow(dt_stem_i) == 0) cli_alert_warning("{.tmp Name} is not in the list for this caclulation.")

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
      cli_abort(c("x" = "There are no such equation type."))
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

      ## adjustment to those in the calculation program (Excel ver)
      ## this is a tentative code. Require "adjust_to_excel" in ellipsis.
      if("adjust_to_excel" %in% names(list(...)) && list(...)["adjust_to_excel"] == TRUE){
        ### 1) Sapporo Todomatsu
        if(Name == "札幌トドマツ" & D >= 90 & D < 92){
          v1 <- calcStemVolumeAdj(dt_stem_i, 90, H, eq_type = 1, adj = "3w")
          v2 <- 0.247459 * H + 0.2424
          V <- linearImpute(90, v1, 92, v2, D)
        }
        if(Name == "札幌トドマツ" & D >= 92 & D < 94){
          wch_replace <- which(D >= 92 & D < 94)
          v1 <- 0.247459 * H + 0.2424
          v2 <- 0.259188 * H + 0.2648
          V  <- linearImpute(92, v1, 94, v2, D)
        }
        if(Name == "札幌トドマツ" & D >= 94 & D < 96){
          wch_replace <- which(D >= 94 & D < 96)
          v1 <- 0.259188 * H + 0.2648
          v2 <- 0.274247 * H + 0.2993
          V  <- linearImpute(94, v1, 96, v2, D)
        }
        if(Name == "札幌トドマツ" & D >= 96 & D < 98){
          wch_replace <- which(D >= 96 & D < 98)
          v1 <- 0.274247 * H + 0.2993
          v2 <- 0.290558 * H + 0.3324
          V  <- linearImpute(96, v1, 98, v2, D)
        }
        if(Name == "札幌トドマツ" & D >= 98 & D < 100){
          wch_replace <- which(D >= 98 & D < 100)
          v1 <- 0.290558 * H + 0.3324
          v2 <- calcStemVolumeAdj(dt_stem_s, 100, H, eq_type = 1, adj = "5w") # use revised coef
          V  <- linearImpute(98, v1, 100, v2, D)
        }
        ### 2) Kochi Tennen-sugi
        # if(Name == "高知天然スギ" & any(D >= 94 & D < 96 & H < 21.5)){ # omitted
        #   wch_replace <- which(D >= 94 & D < 96 & H < 21.5)
        #   v1 <- calcStemVolumeAdj(dt_stem_i, 94, H, eq_type = 1, adj = "3w")
        #   v2 <- 0.253714 * H + 0.1063
        #   V  <- linearImpute(94, v1, 96, v2, D)
        # }
        if(Name == "高知天然スギ" & D >= 96 & D < 98 & H < 21.5){
          wch_replace <- which(D >= 96 & D < 98 & H < 21.5)
          v1 <- 0.253714 * H + 0.1063
          v2 <- calcStemVolumeAdj(dt_stem_i, 98, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(96, v1, 98, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 96 & D < 98 & H >= 29.5){
          wch_replace <- which(D >= 96 & D < 98 & H >= 29.5)
          v1 <- calcStemVolumeAdj(dt_stem_i, 96, H, eq_type = 1, adj = "5w")
          v2 <- 0.249273 * H + 0.4434
          V  <- linearImpute(96, v1, 98, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 98 & D < 100 & H >= 29.5){
          wch_replace <- which(D >= 98 & D < 100 & H >= 29.5)
          v1 <- 0.249273 * H + 0.4434
          v2 <- 0.269714 * H + 0.1319
          V  <- linearImpute(98, v1, 100, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 100 & D < 102 & H >= 29.5){
          wch_replace <- which(D >= 100 & D < 102 & H >= 29.5)
          v1 <- 0.269714 * H + 0.1319
          v2 <- calcStemVolumeAdj(dt_stem_i, 102, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(100, v1, 102, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 102 & D < 104){
          wch_replace <- which(D >= 102 & D < 104)
          v1 <- calcStemVolumeAdj(dt_stem_i, 102, H, eq_type = 1, adj = "5w")
          v2 <- 0.0007 * H^2 + 0.2739 * H - 0.0732
          V  <- linearImpute(102, v1, 104, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 104 & D < 106 & H < 22.5){
          wch_replace <- which(D >= 104 & D < 106 & H < 22.5)
          v1 <- 0.0007 * H^2 + 0.2739 * H - 0.0732
          v2 <- 0.314 * H - 0.467
          V  <- linearImpute(104, v1, 106, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 104 & D < 106 & H >= 22.5 & H < 39.5){
          wch_replace <- which(D >= 104 & D < 106 & H >= 22.5 & H < 39.5)
          v1 <- 0.0007 * H^2 + 0.2739 * H - 0.0732
          v2 <- calcStemVolumeAdj(dt_stem_i, 106, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(104, v1, 106, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 104 & D < 106 & H >= 39.5){
          wch_replace <- which(D >= 104 & D < 106 & H >= 39.5)
          v1 <- 0.0007 * H^2 + 0.2739 * H - 0.0732
          v2 <- 0.354 * H - 1.7846
          V  <- linearImpute(104, v1, 106, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 106 & D < 108 & H < 22.5){
          wch_replace <- which(D >= 106 & D < 108 & H < 22.5)
          v1 <- 0.314 * H - 0.467
          v2 <- 0.327 * H - 0.636
          V  <- linearImpute(106, v1, 108, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 106 & D < 108 &  H >= 22.5 & H < 39.5){
          wch_replace <- which(D >= 106 & D < 108 &  H >= 22.5 & H < 39.5)
          v1 <- calcStemVolumeAdj(dt_stem_i, 106, H, eq_type = 1, adj = "5w")
          v2 <- calcStemVolumeAdj(dt_stem_i, 108, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(106, v1, 108, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 106 & D < 108 & H >= 39.5){
          wch_replace <- which(D >= 106 & D < 108 & H >= 39.5)
          v1 <- 0.354 * H - 1.7846
          v2 <- calcStemVolumeAdj(dt_stem_i, 108, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(106, v1, 108, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 108 & D < 110 & H < 22.5){
          wch_replace <- which(D >= 108 & D < 110 & H < 22.5)
          v1 <- 0.327 * H - 0.636
          v2 <- calcStemVolumeAdj(dt_stem_i, 110, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(108, v1, 110, v2, D)
        }
        if(Name == "高知天然スギ" & D >= 108 & D < 110 & H >= 22.5){
          wch_replace <- which(D >= 108 & D < 110 & H >= 22.5)
          v1 <- calcStemVolumeAdj(dt_stem_i, 108, H, eq_type = 1, adj = "5w")
          v2 <- calcStemVolumeAdj(dt_stem_i, 110, H, eq_type = 1, adj = "5w")
          V  <- linearImpute(108, v1, 110, v2, D)
        }
        ### 3) Aomori Akamatsu
        if(Name == "青森アカマツ" & D > 38 & D <= 46){
          wch_replace <- which(D > 38 & D <= 46)
          v1 <- calcStemVolumeAdj(dt_stem_i, 38, H, eq_type = 1, adj = "None")
          v2 <- calcStemVolumeAdj(dt_stem_i, 46, H, eq_type = 1, adj = "3w")
          V  <- linearImpute(38, v1, 46, v2, D)
        }
        if(Name == "青森アカマツ" & D > 46 & D < 56){
          wch_replace <- which(D > 46 & D < 56)
          v1 <- calcStemVolumeAdj(dt_stem_i, 46, H, eq_type = 1, adj = "3w")
          v2 <- calcStemVolumeAdj(dt_stem_i, 56, H, eq_type = 1, adj = "None")
          V  <- linearImpute(46, v1, 56, v2, D)
        }
        ### 4) Aomori Koyoju, Kochi Koyoju
        if(Name == "青森広葉樹"){
          if(D < 68){
            V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "None")
          } else {
            V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 1, adj = "5w")
          }
        }
        if(Name == "高知広葉樹"){
          if(D < 58){
            V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 2, adj = "3w")
          } else {
            V <- calcStemVolumeAdj(dt_stem_i, D, H, eq_type = 2, adj = "5w")
          }
        }
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
      f_d <- 0.5 - 0.0008*D + 0.421*exp(-0.12*D) + 10^-10 # to avoid floating-point error
      f_h <- 0.61 - 0.0055*H + 5.48*exp(-1.025*H) + 10^-10 # to avoid floating-point error
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
    if("stop_if_NA" %in% names(list(...)) && list(...)[["stop_if_NA"]] == TRUE){
      cli_abort(c("{.tmp Name} should mach one in a list from {.help [{.fun volumeName}](stemv::volumeName)}.",
                  "x" = "{.tmp2 {.var {Name}}} is not in a list of stem volume equation."))
    } else {
      cli_alert_warning("{.tmp2 {.var {Name}}} is not in a list, use {.fun stemv::volumeName}.")
      V <- NA
    }
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
#' 細田ら (2010) 現行立木幹材積表と材積式による計算値との相違およびその修正方法. 森林計画学会誌 44: 23-39.
#' (Kazuo HOSODA, Yasushi MITSUDA and Toshiro IEHARA 2010:
#' "Differences between the present stem volume tables and the values of the volume equations, and their correction"
#' Jpn. J. For. Plann. 44:23-39.)
#'
#' There are several differences compared with the description in the original article and
#' computational program developed based on the same article as follows:
#'
#' - The significant digits for the parameters of several equations are different
#'  (e.g., `函館エゾマツ`, `東京広葉樹`, `長野カラマツ`).
#' - The additional revisions on `札幌トドマツ`,  `高知天然スギ`, and `青森アカマツ` are not reflected.
#' - The range for 5 points moving average adjustment is different for `青森広葉樹` and `高知広葉樹`.
#' - The error in floating-point numbers can cause the differences in rounding calculation for `北海道針葉樹`.
#'
#' @param Name  Japanese character for identifying the calculation methods (e.g., `東京スギ`).
#'   This is expected from [volumeName()].
#' @param D  Numeric. DBH in cm.
#' @param H  Numeric. Tree height in m.
#' @param ... Additional arguments.
#'
#'  These can be:
#'  * `stop_if_NA`: logical. If TRUE, cause error when NAs/NaNs/Inf/-Inf are provided for `Name`/`D`/`H`;
#'    otherwise return NAs.
#'  * `list_data`: data from [getStemCoefficientsRound()] if preferred.
#'  * `adjust_to_excel`: logical. Calculate volume in the same way as program.
#' @return  Calculated stem volume
#'
#' @references \url{https://doi.org/10.20659/jjfp.44.2_23}
#' @importFrom dplyr %>% filter mutate group_by summarize arrange
#' @importFrom purrr map
#' @importFrom utils tail
#' @importFrom cli cli_abort cli_alert_info cli_progress_bar cli_progress_update cli_div
#' @export
#'
stemVolume <- function(Name, D, H, ...){
  ## check D, H --
  cli_div(theme = list(.tmp = list(color = "yellow4", "font-style" = "italic", "font-weight" = "bold"), # set color
                       .tmp2 = list(color = "blue", "font-style" = "italic")))
  if("stop_if_NA" %in% names(list(...)) && list(...)[["stop_if_NA"]] == TRUE){
    if(any(is.na(Name))) cli_abort(c("NA/NaNs should be removed from {.tmp Name}.",
                                     "x" = "There are {sum(is.na(Name))} NA/NaNs in {.tmp Name}."))
    if(any(is.na(D))) cli_abort(c("NA/NaNs should be removed from {.tmp D}.",
                                  "x" = "There are {sum(is.na(D))} NA/NaNs in {.tmp D}."))
    if(any(is.na(H))) cli_abort(c("NA/NaNs should be removed from {.tmp H}.",
                                  "x" = "There are {sum(is.na(H))} NA/NaNs in {.tmp H}."))
    if(any(is.infinite(D))) cli_abort(c("Inf/-Inf should be removed from {.tmp D}.",
                                  "x" = "There are {sum(is.infinite(D))} Inf/-Inf in {.tmp D}."))
    if(any(is.infinite(H))) cli_abort(c("Inf/-Inf should be removed from {.tmp H}.",
                                        "x" = "There are {sum(is.infinite(H))} Inf/-Inf in {.tmp H}."))
  } else {
    if(any(is.na(Name))) cli_alert_warning("There are {sum(is.na(Name))} NA/NaNs in {.tmp Name}.")
    if(any(is.na(D))) cli_alert_warning("There are {sum(is.na(D))} NA/NaNs in {.tmp D}.")
    if(any(is.na(H))) cli_alert_warning("There are {sum(is.na(H))} NA/NaNs in {.tmp H}.")
    if(any(is.infinite(D))) cli_alert_warning("There are {sum(is.infinite(D))} Inf/-Inf in {.tmp D}.")
    if(any(is.infinite(H))) cli_alert_warning("There are {sum(is.infinite(H))} Inf/-Inf in {.tmp H}.")
  }

  ## check length ---
  if(length(Name) == 1 & length(D) > 1){  # in case single Name is provided.
    cli_alert_info("{.tmp Name} is a single, but {.tmp D} is multiple. {.tmp Name} is recycled to length {length(D)}.")
    Name <- rep(Name, length(D))
  }
  if(!all(length(Name) == c(length(D), length(H)))){
    cli_abort(c("Different data length of Name/D/H, but shoud be the same.",
                "x" = paste0("The length of {.tmp Name} is {length(Name)}, ",
                             "the length of {.tmp D} is {length(D)}, ",
                             "the length of {.tmp H} is {length(H)}.")))
  }

  ## check the class --
  if(!is.numeric(D) & !all(is.na(D))) cli_abort(c("{.tmp D} should be numeric.",
                                                  "x" = "class {.tmp D} is {.cls {class(D)}}."))
  if(!is.numeric(H) & !all(is.na(H))) cli_abort(c("{.tmp H} should be numeric.",
                                                  "x" = "class {.tmp H} is {.cls {class(H)}}."))

  ## set progress bar if data is large ---
  if(length(D) > 10^5)  cli_progress_bar("Calculating data", total = length(unique(Name)))

  ## adjustment to D & H. if either is NA, return NA ---
  if(any(is.infinite(D))) D[is.infinite(D)] <- NA
  if(any(is.infinite(H))) H[is.infinite(H)] <- NA
  if(any(D < 0, na.rm = T)){
    cli_alert_warning("There are {sum(D < 0, na.rm = T)} negative values (< 0) in {.tmp D}, adjusting to 0.")
    D[D < 0] <- 0
  }
  if(any(H < 0, na.rm = T)){
    cli_alert_warning("There are {sum(H < 0, na.rm = T)} negative values (< 0) in {.tmp H}, adjusting to 0.")
    H[H < 0] <- 0
  }


  ## prepare V for all measurements ---
  V <- rep(NA, length(D))

  ## data for coefficients ---
  ### use external coefficients if provided in ellipsis
  ### otherwise call the function to get coefficient
  if("list_data" %in% names(list(...))){
    list_data <- list(...)[["list_data"]]
  } else {
    list_data <- getStemCoefficients()
  }
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
    cli_alert_warning("There are {.tmp Name}s that are not listed in this caclulation.")

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
        cli_abort(c("{.tmp eq_type_i} should be integer between 1 and 4.", "x" = "There are no such equation type."))
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
      if(ls_Name[i] == "上屋久天然スギ" & any(D_i >= 61, na.rm = T)){
        wch_replace <- which(D_i >= 61)
        V_i[wch_replace] = 0.1513 + 0.8898 * V_i[wch_replace]
      }
      if(ls_Name[i] == "上屋久天然スギ" & any(abs(D_i - 61) < 5, na.rm = T)){
        dt_stem_s <- dt_stem_i |> rbind(dt_stem_i[2,]) |> arrange(D_lower)
        dt_stem_s$D_upper[2] <- 61
        dt_stem_s$D_lower[3] <- 61
        wch_replace <- which(abs(D_i - 61) < 5)
        V_i[wch_replace] = calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 5, adj = "5w")
      }
      if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 95, na.rm = T)){
        dt_stem_s <- dt_stem |> filter(name == "札幌トドマツ") |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        wch_replace <- which(D_i >= 95)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "札幌トドマツ" & any(abs(D_i - 95) < 5, na.rm = T)){
        dt_stem_s <- dt_stem |> filter(name == "札幌トドマツ") |> tail(1) |> mutate(D_upper = 95) |>
          rbind(dt_stem |> filter(name == "札幌エゾマツ") |> tail(1) |> mutate(D_lower = 95))
        wch_replace <- which(abs(D_i - 95) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_s, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "青森アカマツ" & any(abs(D_i - 47) < 5, na.rm = T)){
        wch_replace <- which(abs(D_i - 47) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "青森広葉樹" & any(abs(D_i - 71) < 5, na.rm = T)){
        wch_replace <- which(abs(D_i - 71) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "高知天然スギ" & any(abs(D_i - 101) < 5, na.rm = T)){
        wch_replace <- which(abs(D_i - 101) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
      }
      if(ls_Name[i] == "高知広葉樹" & any(abs(D_i - 61) < 5, na.rm = T)){
        wch_replace <- which(abs(D_i - 61) < 5)
        V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 2, adj = "5w")
      }
      if(ls_Name[i] == "東京スギ" & any(H_i < 11.5 & D_i >= 36 & D_i <= 70, na.rm = T)){
        wch_replace <- which(H_i < 11.5 & D_i >= 36 & D_i <= 70)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 70, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 11.5 & H_i < 12.5 & D_i >= 36 & D_i <= 68, na.rm = T)){
        wch_replace <- which(H_i >= 11.5 & H_i < 12.5 & D_i >= 36 & D_i <= 68)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 68, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 12.5 & H_i < 14.5 & D_i >= 36 & D_i <= 62, na.rm = T)){
        wch_replace <- which(H_i >= 12.5 & H_i < 14.5 & D_i >= 36 & D_i <= 62)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 62, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 14.5 & H_i < 16.5 & D_i >= 36 & D_i <= 60, na.rm = T)){
        wch_replace <- which(H_i >= 14.5 & H_i < 16.5 & D_i >= 36 & D_i <= 60)
        V_i[wch_replace] <- linearImputeVolume(dt_stem_i, 36, 60, D_i[wch_replace], H_i[wch_replace], eq_type = 1)
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 16.5 & H_i < 18.5 & D_i >= 40 & D_i <= 60, na.rm = T)){
        wch_replace <- which(H_i >= 16.5 & H_i < 18.5 & D_i >= 40 & D_i <= 60)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 60, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 60, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 18.5 & H_i < 19.5 & D_i >= 40 & D_i <= 58, na.rm = T)){
        wch_replace <- which(H_i >= 18.5 & H_i < 19.5 & D_i >= 40 & D_i <= 58)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 58, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 58, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 19.5 & H_i < 20.5 & D_i >= 40 & D_i <= 54, na.rm = T)){
        wch_replace <- which(H_i >= 19.5 & H_i < 20.5 & D_i >= 40 & D_i <= 54)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 54, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 54, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 20.5 & H_i < 21.5 & D_i >= 40 & D_i <= 52, na.rm = T)){
        wch_replace <- which(H_i >= 20.5 & H_i < 21.5 & D_i >= 40 & D_i <= 52)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 52, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 52, v2, D_i[wch_replace])
      }
      if(ls_Name[i] == "東京スギ" & any(H_i >= 21.5 & H_i < 22.5 & D_i >= 40 & D_i <= 50, na.rm = T)){
        wch_replace <- which(H_i >= 21.5 & H_i < 22.5 & D_i >= 40 & D_i <= 50)
        v1 <- calcStemVolumeMulti(dt_stem_i, 40, H_i[wch_replace], eq_type = 1, adj = "3w")
        v2 <- calcStemVolumeMulti(dt_stem_i, 50, H_i[wch_replace], eq_type = 1, adj = "None")
        V_i[wch_replace]  <- linearImpute(40, v1, 50, v2, D_i[wch_replace])
      }

      ## adjustment to those in the calculation program (Excel ver)
      ## this is a tentative code. Require "adjust_to_excel" in ellipsis.
      if("adjust_to_excel" %in% names(list(...)) && list(...)[["adjust_to_excel"]] == TRUE){
        ### 1) Sapporo Todomatsu
        if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 90 & D_i < 92, na.rm = T)){
          wch_replace <- which(D_i >= 90 & D_i < 92)
          v1 <- calcStemVolumeMulti(dt_stem_i, 90, H_i[wch_replace], eq_type = 1, adj = "3w")
          v2 <- 0.247459 * H_i[wch_replace] + 0.2424
          V_i[wch_replace]  <- linearImpute(90, v1, 92, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 92 & D_i < 94, na.rm = T)){
          wch_replace <- which(D_i >= 92 & D_i < 94)
          v1 <- 0.247459 * H_i[wch_replace] + 0.2424
          v2 <- 0.259188 * H_i[wch_replace] + 0.2648
          V_i[wch_replace]  <- linearImpute(92, v1, 94, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 94 & D_i < 96, na.rm = T)){
          wch_replace <- which(D_i >= 94 & D_i < 96)
          v1 <- 0.259188 * H_i[wch_replace] + 0.2648
          v2 <- 0.274247 * H_i[wch_replace] + 0.2993
          V_i[wch_replace]  <- linearImpute(94, v1, 96, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 96 & D_i < 98, na.rm = T)){
          wch_replace <- which(D_i >= 96 & D_i < 98)
          v1 <- 0.274247 * H_i[wch_replace] + 0.2993
          v2 <- 0.290558 * H_i[wch_replace] + 0.3324
          V_i[wch_replace]  <- linearImpute(96, v1, 98, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "札幌トドマツ" & any(D_i >= 98 & D_i < 100, na.rm = T)){
          wch_replace <- which(D_i >= 98 & D_i < 100)
          v1 <- 0.290558 * H_i[wch_replace] + 0.3324
          v2 <- calcStemVolumeMulti(dt_stem_s, 100, H_i[wch_replace], eq_type = 1, adj = "5w") # use revised coef
          V_i[wch_replace]  <- linearImpute(98, v1, 100, v2, D_i[wch_replace])
        }
        ### 2) Kochi Tennen-sugi
        # if(ls_Name[i] == "高知天然スギ" & any(D_i >= 94 & D_i < 96 & H_i < 21.5)){ # omitted
        #   wch_replace <- which(D_i >= 94 & D_i < 96 & H_i < 21.5)
        #   v1 <- calcStemVolumeMulti(dt_stem_i, 94, H_i[wch_replace], eq_type = 1, adj = "3w")
        #   v2 <- 0.253714 * H_i[wch_replace] + 0.1063
        #   V_i[wch_replace]  <- linearImpute(94, v1, 96, v2, D_i[wch_replace])
        # }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 96 & D_i < 98 & H_i < 21.5, na.rm = T)){
          wch_replace <- which(D_i >= 96 & D_i < 98 & H_i < 21.5)
          v1 <- 0.253714 * H_i[wch_replace] + 0.1063
          v2 <- calcStemVolumeMulti(dt_stem_i, 98, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(96, v1, 98, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 96 & D_i < 98 & H_i >= 29.5, na.rm = T)){
          wch_replace <- which(D_i >= 96 & D_i < 98 & H_i >= 29.5)
          v1 <- calcStemVolumeMulti(dt_stem_i, 96, H_i[wch_replace], eq_type = 1, adj = "5w")
          v2 <- 0.249273 * H_i[wch_replace] + 0.4434
          V_i[wch_replace]  <- linearImpute(96, v1, 98, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 98 & D_i < 100 & H_i >= 29.5, na.rm = T)){
          wch_replace <- which(D_i >= 98 & D_i < 100 & H_i >= 29.5)
          v1 <- 0.249273 * H_i[wch_replace] + 0.4434
          v2 <- 0.269714 * H_i[wch_replace] + 0.1319
          V_i[wch_replace]  <- linearImpute(98, v1, 100, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 100 & D_i < 102 & H_i >= 29.5, na.rm = T)){
          wch_replace <- which(D_i >= 100 & D_i < 102 & H_i >= 29.5)
          v1 <- 0.269714 * H_i[wch_replace] + 0.1319
          v2 <- calcStemVolumeMulti(dt_stem_i, 102, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(100, v1, 102, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 102 & D_i < 104, na.rm = T)){
          wch_replace <- which(D_i >= 102 & D_i < 104)
          v1 <- calcStemVolumeMulti(dt_stem_i, 102, H_i[wch_replace], eq_type = 1, adj = "5w")
          v2 <- 0.0007 * H_i[wch_replace]^2 + 0.2739 * H_i[wch_replace] - 0.0732
          V_i[wch_replace]  <- linearImpute(102, v1, 104, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 104 & D_i < 106 & H_i < 22.5, na.rm = T)){
          wch_replace <- which(D_i >= 104 & D_i < 106 & H_i < 22.5)
          v1 <- 0.0007 * H_i[wch_replace]^2 + 0.2739 * H_i[wch_replace] - 0.0732
          v2 <- 0.314 * H_i[wch_replace] - 0.467
          V_i[wch_replace]  <- linearImpute(104, v1, 106, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 104 & D_i < 106 & H_i >= 22.5 & H_i < 39.5, na.rm = T)){
          wch_replace <- which(D_i >= 104 & D_i < 106 & H_i >= 22.5 & H_i < 39.5)
          v1 <- 0.0007 * H_i[wch_replace]^2 + 0.2739 * H_i[wch_replace] - 0.0732
          v2 <- calcStemVolumeMulti(dt_stem_i, 106, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(104, v1, 106, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 104 & D_i < 106 & H_i >= 39.5, na.rm = T)){
          wch_replace <- which(D_i >= 104 & D_i < 106 & H_i >= 39.5)
          v1 <- 0.0007 * H_i[wch_replace]^2 + 0.2739 * H_i[wch_replace] - 0.0732
          v2 <- 0.354 * H_i[wch_replace] - 1.7846
          V_i[wch_replace]  <- linearImpute(104, v1, 106, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 106 & D_i < 108 & H_i < 22.5, na.rm = T)){
          wch_replace <- which(D_i >= 106 & D_i < 108 & H_i < 22.5)
          v1 <- 0.314 * H_i[wch_replace] - 0.467
          v2 <- 0.327 * H_i[wch_replace] - 0.636
          V_i[wch_replace]  <- linearImpute(106, v1, 108, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 106 & D_i < 108 &  H_i >= 22.5 & H_i < 39.5, na.rm = T)){
          wch_replace <- which(D_i >= 106 & D_i < 108 &  H_i >= 22.5 & H_i < 39.5)
          v1 <- calcStemVolumeMulti(dt_stem_i, 106, H_i[wch_replace], eq_type = 1, adj = "5w")
          v2 <- calcStemVolumeMulti(dt_stem_i, 108, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(106, v1, 108, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 106 & D_i < 108 & H_i >= 39.5, na.rm = T)){
          wch_replace <- which(D_i >= 106 & D_i < 108 & H_i >= 39.5)
          v1 <- 0.354 * H_i[wch_replace] - 1.7846
          v2 <- calcStemVolumeMulti(dt_stem_i, 108, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(106, v1, 108, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 108 & D_i < 110 & H_i < 22.5, na.rm = T)){
          wch_replace <- which(D_i >= 108 & D_i < 110 & H_i < 22.5)
          v1 <- 0.327 * H_i[wch_replace] - 0.636
          v2 <- calcStemVolumeMulti(dt_stem_i, 110, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(108, v1, 110, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "高知天然スギ" & any(D_i >= 108 & D_i < 110 & H_i >= 22.5, na.rm = T)){
          wch_replace <- which(D_i >= 108 & D_i < 110 & H_i >= 22.5)
          v1 <- calcStemVolumeMulti(dt_stem_i, 108, H_i[wch_replace], eq_type = 1, adj = "5w")
          v2 <- calcStemVolumeMulti(dt_stem_i, 110, H_i[wch_replace], eq_type = 1, adj = "5w")
          V_i[wch_replace]  <- linearImpute(108, v1, 110, v2, D_i[wch_replace])
        }
        ### 3) Aomori Akamatsu
        if(ls_Name[i] == "青森アカマツ" & any(D_i > 38 & D_i <= 46, na.rm = T)){
          wch_replace <- which(D_i > 38 & D_i <= 46)
          v1 <- calcStemVolumeMulti(dt_stem_i, 38, H_i[wch_replace], eq_type = 1, adj = "None")
          v2 <- calcStemVolumeMulti(dt_stem_i, 46, H_i[wch_replace], eq_type = 1, adj = "3w")
          V_i[wch_replace]  <- linearImpute(38, v1, 46, v2, D_i[wch_replace])
        }
        if(ls_Name[i] == "青森アカマツ" & any(D_i > 46 & D_i < 56, na.rm = T)){
          wch_replace <- which(D_i > 46 & D_i < 56)
          v1 <- calcStemVolumeMulti(dt_stem_i, 46, H_i[wch_replace], eq_type = 1, adj = "3w")
          v2 <- calcStemVolumeMulti(dt_stem_i, 56, H_i[wch_replace], eq_type = 1, adj = "None")
          V_i[wch_replace]  <- linearImpute(46, v1, 56, v2, D_i[wch_replace])
        }
        ### 4) Aomori Koyoju, Kochi Koyoju
        if(ls_Name[i] == "青森広葉樹"){
          wch_replace <- which(D_i < 68)
          V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "None")
          wch_replace <- which(D_i >= 68)
          V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 1, adj = "5w")
        }
        if(ls_Name[i] == "高知広葉樹"){
          wch_replace <- which(D_i < 58)
          V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 2, adj = "3w")
          wch_replace <- which(D_i >= 58)
          V_i[wch_replace] <- calcStemVolumeMulti(dt_stem_i, D_i[wch_replace], H_i[wch_replace], eq_type = 2, adj = "5w")
        }
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
        f_d <- 0.5 - 0.0008*D_i+ 0.421*exp(-0.12*D_i) + 10^-10 # to avoid floating-point error
        f_h <- 0.61 - 0.0055*H_i + 5.48*exp(-1.025*H_i) + 10^-10 # to avoid floating-point error
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

      ## does not much any stem volume equation ---
    } else {
      if("stop_if_NA" %in% names(list(...)) && list(...)[["stop_if_NA"]] == TRUE){
        cli_abort(c("{.tmp Name} should mach one in a list from {.help [{.fun volumeName}](stemv::volumeName)}.",
                    "x" = "{.tmp2 {.var {ls_Name[i]}}} is not in a list of stem volume equation."))
      } else {
        cli_alert_warning("{.tmp2 {.var {ls_Name[i]}}} is not in a list, use {.fun stemv::volumeName}.")
      }
    }

    ## combine the calculated volume ---
    V[wch_i] <- V_i

    if(length(D) > 10^5) cli_progress_update()
  }
  if(any(is.na(V))) cli_alert_warning("The calculation of {.tmp V} contatins {sum(is.na(V))} NA/NaNs.")
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
#' @param ...  Additional arguments. `name_invalid` and `stop_if_NA` are supported.
#'
#'  These can be:
#'  * `stop_if_NA`: logical. If TRUE, cause error when NAs/NaNs/NULL or other invalid region name is provided for `Region`.
#'  * `name_invalid`: character (or NA/NULL), which is used as return when invalid region name is provided for `Region`.
#'   This is not used if `stop_if_NA` is `TRUE`.
#' @return  Name for identifying and calculating stem volume in [stemVolume()]
#'
#' @importFrom dplyr %>% filter mutate group_by summarize
#' @importFrom stringr str_replace str_detect
#' @importFrom stringi stri_trans_general
#' @importFrom cli cli_abort cli_alert_warning cli_div
#' @export

volumeNameSingle <- function(Region, Spp, RS = NULL, ...){
  # determine returns for NAs & invalid names provided
  if("name_invalid" %in% names(list(...))){
    name_invalid <- list(...)[["name_invalid"]]
  } else {
    name_invalid <- NA  # return NA by default
  }
  if("stop_if_NA" %in% names(list(...))){
    stop_if_NA <- list(...)[["stop_if_NA"]]
  } else {
    stop_if_NA <- FALSE  # do not cause error if NA by default
  }

  ### check length ---
  cli_div(theme = list(.tmp = list(color = "yellow4", "font-style" = "italic", "font-weight" = "bold"), # set color
                       .tmp2 = list(color = "blue", "font-style" = "italic")))
  if(any(c(length(Region) != 1, length(Region) != 1)))
    cli_abort(c("x" = "Only single observation is allowed for Region/Spp."))
  ## stop_if_NA = TRUE & Spp = NA, then cause error
  if(stop_if_NA){
    if(is.na(Spp)) cli_abort(c("{.tmp Spp} should be character.", "x" = "{.tmp Spp} is NA/NaN."))
    if(!is.character(Spp)) cli_abort(c("{.tmp Spp} should be character.", "x" = "class {.tmp Spp} is {.cls {class(Spp)}}."))
  }
  ## stop_if_NA = FALSE & Spp = NA, then return NA
  if(!stop_if_NA){
    if(is.na(Spp)){
      cli_alert_warning(c("{.tmp Spp} should be character, return {.tmp2 {.var {name_invalid}}}."))
      return(name_invalid)
    }
    if(!is.character(Spp)) cli_abort(c("{.tmp Spp} should be character.", "x" = "class {.tmp Spp} is {.cls {class(Spp)}}."))
  }

  # prepare for Region ---
  name_region <- Region
  name_region <- str_replace(name_region, "(府|県)$", "")
  name_region <- str_replace(name_region, "東京都", "東京") # to avoid "京都" to "都"
  # prepare for Spp ---
  name_spp <- Spp
  if(str_detect(name_spp, pattern = "\\p{Katakana}|\\p{Han}"))
    name_spp <- stri_trans_general(name_spp, "Halfwidth-Fullwidth")
  if(str_detect(name_spp, pattern = "\\p{Hiragana}"))
    name_spp <- stri_trans_general(name_spp, "Hiragana-Katakana")
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
    } else if(name_spp %in% c(RS$s_x_Sugi, RS$s_x_TenSugi)){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      cli_abort("No such species in this region (Asahikawa).")
    }
    ## 2. Kitami
  } else if (name_region %in% RS$r_2_Kitami){
    if(name_spp %in% RS$s_2_Todomatsu){
      name_RS <- "北見トドマツ"
    } else if(name_spp %in% RS$s_2_Ezomatsu){
      name_RS <- "北見エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% c(RS$s_x_Sugi, RS$s_x_TenSugi)){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      cli_abort("No such species in this region (Kitami).")
    }
    ## 3. Obihiro
  } else if (name_region %in% RS$r_3_Obihiro){
    if(name_spp %in% RS$s_3_Todomatsu){
      name_RS <- "帯広トドマツ"
    } else if(name_spp %in% RS$s_3_Ezomatsu){
      name_RS <- "帯広エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% c(RS$s_x_Sugi, RS$s_x_TenSugi)){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      cli_abort("No such species in this region (Obihiro).")
    }
    ## 4. Sapporo
  } else if (name_region %in% RS$r_4_Sapporo){
    if(name_spp %in% RS$s_4_Todomatsu){
      name_RS <- "札幌トドマツ"
    } else if(name_spp %in% RS$s_4_Ezomatsu){
      name_RS <- "札幌エゾマツ"
    } else if(name_spp %in% RS$s_x_Karamatsu){
      name_RS <- "北海道カラマツ"
    } else if(name_spp %in% c(RS$s_x_Sugi, RS$s_x_TenSugi)){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      cli_abort("No such species in this region (Sapporo).")
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
    } else if(name_spp %in% c(RS$s_x_Sugi, RS$s_x_TenSugi)){
      name_RS <- "青森スギ"
    } else if(any(str_detect(name_spp, RS$s_x_ConiferOther))){
      name_RS <- "北海道針葉樹"
    } else if(any(str_detect(name_spp, RS$s_x_BroadleafOther))){
      name_RS <- "北海道広葉樹"
    } else {
      cli_abort("No such species in this region (Hakodate).")
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
      cli_abort("No such species in this region (Aomori).")
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
      cli_abort("No such species in this region (Akita).")
    }
    ## 8. Maebashi
  } else if (name_region %in% RS$r_8_Maebashi){
    if(name_spp %in% RS$s_8_Sugi & name_region %in% c("会津新潟", "新潟", "会津地方")){
      name_RS <- "会津新潟スギ"
    } else if(name_spp %in% RS$s_8_Sugi){
      name_RS <- "前橋スギ"
    } else if(name_spp %in% RS$s_8_Hinoki){
      name_RS <- "前橋ヒノキ"
    } else if(name_spp %in% RS$s_8_Karamatsu){
      name_RS <- "前橋カラマツ"
    } else if(name_spp %in% RS$s_8_Akamatsu & name_region %in% c("会津新潟", "新潟", "会津地方")){
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
      cli_abort("No such species in this region (Maebashi).")
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
      cli_abort("No such species in this region (Tokyo).")
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
      cli_abort("No such species in this region (Nagano).")
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
      cli_abort("No such species in this region (Nagoya).")
    }
    ## 12. Osaka
  } else if (name_region %in% RS$r_12_Osaka){
    if(name_spp %in% RS$s_12_Sugi){
      name_RS <- "大阪スギ"
    } else if(name_spp %in% RS$s_12_TenSugi & !(name_region %in% c("山陰", "鳥取", "島根", "石川", "福井",
                                                                     "中丹地方", "丹後地方", "但馬地方"))){
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
      cli_abort("No such species in this region (Osaka).")
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
      cli_abort("No such species in this region (Kochi).")
    }
    ## 14. Kumamoto
  } else if (name_region %in% RS$r_14_Kumamoto){
    if(name_spp %in% c(RS$s_14_Sugi, RS$s_14_TenSugi) & name_region %in% c("飫肥", "飫肥地方")){
      name_RS <- "飫肥スギ"
    } else if(name_spp %in% RS$s_14_Sugi){
      name_RS <- "熊本スギ"
    } else if(name_spp %in% RS$s_14_TenSugi & !(name_region %in% c("下屋久", "下屋久地方", "上屋久", "上屋久地方"))){
      name_RS <- "熊本スギ"
    } else if(name_spp %in% RS$s_14_TenSugi & name_region %in% c("下屋久", "下屋久地方")){
      name_RS <- "下屋久天然スギ"
    } else if(name_spp %in% RS$s_14_TenSugi){
      name_RS <- "上屋久天然スギ"
    } else if(name_spp %in% RS$s_14_Hinoki){
      name_RS <- "熊本ヒノキ"
    } else if(name_spp %in% RS$s_14_Akamatsu & name_region %in% c("霧島", "霧島地方")){
      name_RS <- "霧島アカマツ"
    } else if(name_spp %in% RS$s_14_TenAkamatsu & name_region %in% c("霧島", "霧島地方")){
      name_RS <- "霧島天然アカマツ"
    } else if(name_spp %in% c(RS$s_14_Akamatsu, RS$s_14_TenAkamatsu)){
      name_RS <- "熊本アカマツ"
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
      cli_abort("No such species in this region (Kumamoto).")
    }
  } else {
    if(stop_if_NA){
      cli_abort(c("x" = "{.tmp2 {.var {Region}}} is not in the region list. There are no such region."))
    } else {
      cli_alert_warning(c("{.tmp2 {.var {Region}}} is not in the region list, return {.tmp2 {.var {name_invalid}}}."))
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
#' @param ...  Additional arguments. `name_invalid` and `stop_if_NA` are supported.
#'
#'  These can be:
#'  * `stop_if_NA`: logical. If TRUE, cause error when NAs/NaNs/NULL or other invalid region name is provided for `Region`.
#'  * `name_invalid`: character (or NA/NULL), which is used as return when invalid region name is provided for `Region`.
#'   This is not used if `stop_if_NA` is `TRUE`.
#' @return Name for identifying and calculating stem volume in [stemVolume()]
#'
#' @importFrom purrr pmap
#' @importFrom dplyr %>% filter mutate group_by summarize
#' @importFrom cli cli_abort cli_alert_warning cli_div
#' @export
#'
volumeName <- function(Region, Spp, RS = NULL, ...){
  # for single value ---
  if(length(Region) == 1 & length(Spp) == 1){
    Name <- volumeNameSingle(Region, Spp, RS, ...)

  # for multiple values ---
  } else {
    ## check length ---
    cli_div(theme = list(.tmp = list(color = "yellow4", "font-style" = "italic", "font-weight" = "bold"), # set color
                         .tmp2 = list(color = "blue", "font-style" = "italic")))
    if(length(Region) != length(Spp))
      cli_abort(c("The length should be the same for Region/Spp.",
                  "x" = "The length of {.tmp Region} is {length(Region)}, but the length of {.tmp Spp} is {length(Spp)}."))
    ## check class ---
    if(any(!is.na(Spp) & !is.character(Spp)))
      cli_abort(c("{.tmp Spp} should be character.",
                  "x" = "There are {sum(!is.na(Spp) & !is.character(Spp))} {.tmp Spp} that are not {.cls character}."))

    ## unique combinations of Region & Spp ---
    dt_comb <- data.frame(Region = Region, Spp = Spp) |>
      group_by(Region, Spp) |>
      summarize(.groups = "drop")
    ## get Name for each unique comb. ---
    Name_comb <- pmap(dt_comb, volumeNameSingle, RS, ...)
    ## allocate to each Region & Spp ---
    Name <- factor(paste0(Region, Spp),
                   levels = paste0(dt_comb$Region, dt_comb$Spp),
                   labels = unlist(Name_comb)) |>
      as.character()
  }
  return(Name)
}


