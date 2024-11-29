#' # Compute Zone of Influence for points and annotate them to points using SQL
#'
#' Need to check functions (gaussian, exp, linear) and check if this is correct
#' Need to implement argument zoi_limit for exp and gaussian
#'
#' @param input_points `[character]` \cr Name of the input table of points to be annotated (within the connection `con`,
#' a database), in the format `table_name` or `schema_name.table_name`.
#' @param infrastructure_layer `[character]` \cr Name of the infrastructure covariate table for which
#' the zone of influence will be computed (within the connection `con`, a database),
#' in the format `table_name` or `schema_name.table_name`.
#' @param radius `[numeric=100]` \cr Radius or scale of zone of influence, used to calculate the
#' cumulative ZOI and density. The radius represent the distance at which the ZOI vanishes or
#' goes below a given minimum limit value `zoi_limit`. ### CHECK THAT, so far only bartlett
#' @param input_geom `[character]` \cr Name of the geometry column from the `input_points` table.
#' @param infra_geom `[character]` \cr Name of the geometry column from the `infrastructure_layer` table.
#' @param input_id `[character]` \cr Name of a ID column from the `input_points` table.
#' @param output_type `[character]` \cr
#' @param output_column_name `[character]` \cr
#'
#' @export

# declaring the function
calc_zoi_sql <- function(input_points,
                         infrastructure_layer,
                         radius = 100,
                         type = c("circle", "Gauss", "exp_decay", "bartlett", "threshold")[1],
                         zoi_metric = c("cumulative", "nearest")[1],
                         input_id = "id",
                         input_geom = "geom",
                         infra_geom = "geom",
                         output_type = c("cumulative_zoi", "density")[1], # only cumulative_zoi working
                         output_column_name = paste0("infrastructure_layer_", output_type, "_bartlett", radius),
                         limit = 1000000000000000,
                         verbose = FALSE) {

  if(type %in% c("circle", "threshold")) {
    dd_func <- "coalesce(sign(?radius - ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)"
    # radius_div <- 1
  } else {
    if(type %in% c("bartlett", "linear")) {
      dd_func <- "coalesce((1 - ST_Distance(pts.?input_geom, infra.?infra_geom)/?radius), 0)"
      # radius_div <- radius#ifelse(output_type == "cumulative_zoi", 1, radius)
    } else {
      if(type %in% c("exp_decay", "exponential")) {
        dd_func <- "coalesce(exp(- log(1/0.05)/?radius * ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)"
        # dd_func <- "coalesce(exp(- log(1/?zoi_limit)/?radius * ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)"
        # radius_div <- ifelse(output_type == "cumulative_zoi", 1, radius)
      }
    }
  }

  base_query <- paste0("
SELECT pts.?input_points_id,
  ?summary(", dd_func, ") AS ?out_name
  FROM ?input_pts AS pts
    LEFT JOIN ?infra_layer AS infra
    ON ST_DWithin(pts.?input_geom, infra.?infra_geom, ?radius) AND date_part('year', pts.acquisition_time) >= infra.\"StartActiv\"
  GROUP BY pts.?input_points_id
limit ?lim;
")

  qq <- DBI::sqlInterpolate(con, base_query,
                            # dist_decay_func = DBI::SQL(),
                            input_pts = DBI::SQL(input_points),
                            infra_layer = DBI::SQL(infrastructure_layer),
                            radius = DBI::SQL(radius),
                            # radius_div = DBI::SQL(radius_div),
                            input_geom = DBI::SQL(input_geom),
                            infra_geom = DBI::SQL(infra_geom),
                            input_points_id = DBI::SQL(input_id),
                            out_name = DBI::SQL(output_column_name),
                            summary = DBI::SQL(ifelse(zoi_metric == "cumulative", "sum", "max")),
                            lim = limit)
  qq
  if(verbose) print(qq)

  DBI::dbGetQuery(con, qq)
}

# using the function
# source("~/.pgpass")
#
# NinaR::postgreSQLConnect(
#   host = "gisdata-db.nina.no",
#   dbname = "gisdata",
#   username = pg_username,
#   password = pg_password
# )
#
# calc_zoi_cumulative_sql(input_points = "sam_trein_ancillary.use_ava_data_trein_nose_traj26h_oneimpact_env_f",
#                         infrastructure_layer = "sam_env.wind_turbines_no",
#                         radius = 5000,
#                         type = "exp_decay",
#                         input_geom = "pt_geom_end_e33",
#                         infra_geom = "geom",
#                         input_id = "use_ava_data_animals_id",
#                         output_type = "cumulative_zoi",
#                         zoi_metric = "cumulative",
#                         limit = 100)



