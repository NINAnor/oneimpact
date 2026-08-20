#' Compute point-level Zone of Influence (ZOI) values with SQL
#'
#' Computes ZOI-based covariates from an infrastructure layer and annotates them
#' to input points using SQL in a database connection (for example PostgreSQL).
#' Unlike raster-based ZOI functions, this computes values only at input point locations.
#' For optimal use, both `input_points` and `infrastructure_layer` should be spatially indexed,
#' so only nearby points are evaluated and the query is completed fast.
#'
#' So far, only linear decay/bartlett is implemented.
#'
#' @param con `[DBIConnection]` \cr Open database connection used to run the SQL query.
#' @param input_points `[character(1)]` \cr Name of the input points table, in the format
#' `table_name` or `schema_name.table_name`. The table needs to be within
#' the connection `con`, a database.
#' @param infrastructure_layer `[character(1)]` \cr Name of the infrastructure table used to
#' compute ZOI, in the format `table_name` or `schema_name.table_name`.
#' @param radius `[numeric]` \cr Radius used in the ZOI function. Can be a vector
#' to compute multiple output columns.
#' @param type `[character(1)="circle"]{"circle","Gauss","exp_decay","bartlett","threshold","linear","exponential"}` \cr
#' Distance-decay function used to compute ZOI values.
#' @param zoi_metric `[character(1)="cumulative"]{"cumulative","nearest"}` \cr Aggregation metric:
#' `"cumulative"` uses `sum`, and `"nearest"` uses `max`.
#' @param input_id `[character(1)="id"]` \cr Identifier column in `input_points`.
#' @param input_geom `[character(1)="geom"]` \cr Geometry column in `input_points`.
#' @param infra_geom `[character(1)="geom"]` \cr Geometry column in `infrastructure_layer`.
#' @param output_type `[character(1)="cumulative_zoi"]{"cumulative_zoi","density"}` \cr Output
#' interpretation used in naming and downstream use.
#' @param zoi_limit `[numeric(1)=0.05]` \cr Limit value used in decay-based functions (for example
#' exponential decay calibration). Not yet implemented, to be updated.
#' @param output_column_name `[character]` \cr Output column name(s) for computed ZOI values.
#' Length should match `radius` when multiple radii are provided.
#' @param output_table `[character(1)=NULL]` \cr Optional output table name. If provided, the query
#' is prefixed with `CREATE TABLE IF NOT EXISTS ... AS`.
#' @param condition `[character(1)=""]` \cr Optional extra SQL condition appended to the spatial join
#' clause (for example temporal filtering).
#' @param limit `[numeric(1)=1000000000000000]` \cr Row limit applied to the final query.
#' Useful for limiting the number of rows for testing purposes.
#' @param verbose `[logical(1)=FALSE]` \cr Whether to print the generated SQL query.
#'
#' @return A `data.frame` with one row per input point and computed ZOI column(s). The output
#' includes the input point ID column and one or more derived ZOI columns named by
#' `output_column_name`. If `output_table` is provided, table creation is attempted before query
#' execution, and returned result content may depend on database backend behavior.
#'
#' @example examples/calc_zoi_sql_example.R
#'
#' @export
calc_zoi_sql <- function(con,
                         input_points,
                         infrastructure_layer,
                         radius = 100,
                         type = c("circle", "Gauss", "exp_decay", "bartlett", "threshold")[1],
                         zoi_metric = c("cumulative", "nearest")[1],
                         input_id = "id",
                         input_geom = "geom",
                         infra_geom = "geom",
                         output_type = c("cumulative_zoi", "density")[1], # only cumulative_zoi working
                         zoi_limit = 0.05,
                         output_column_name = paste0(infrastructure_layer, "_", output_type, "_bartlett", radius),
                         output_table = NULL,
                         condition = "", #"AND date_part('year', pts.acquisition_time) > infra.\"StartActiv\""
                         limit = 1000000000000000,
                         verbose = FALSE) {

  func_list <- c()

  for(i in seq_along(radius)) {
    if(type %in% c("circle", "threshold")) {
      dd_func <- DBI::sqlInterpolate(con,
                                     "?summary(coalesce(sign(?radius - ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)) AS ?out_col_name",
                                     summary = DBI::SQL(ifelse(zoi_metric == "cumulative", "sum", "max")),
                                     radius = DBI::SQL(radius[i]),
                                     input_geom = DBI::SQL(input_geom),
                                     infra_geom = DBI::SQL(infra_geom),
                                     out_col_name = DBI::SQL(output_column_name)[i])
      # radius_div <- 1
    } else {
      if(type %in% c("bartlett", "linear")) {
        dd_func <- DBI::sqlInterpolate(con,
                                       "?summary(coalesce((1 - ST_Distance(pts.?input_geom, infra.?infra_geom)/?radius), 0)) AS ?out_col_name",
                                       summary = DBI::SQL(ifelse(zoi_metric == "cumulative", "sum", "max")),
                                       radius = DBI::SQL(radius[i]),
                                       input_geom = DBI::SQL(input_geom),
                                       infra_geom = DBI::SQL(infra_geom),
                                       out_col_name = DBI::SQL(output_column_name)[i])
        # radius_div <- radius#ifelse(output_type == "cumulative_zoi", 1, radius)
      } else {
        if(type %in% c("exp_decay", "exponential")) {
          dd_func <- DBI::sqlInterpolate(con,
                                         "?summary(coalesce(exp(- log(1/?zoi_limit)/?radius * ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)) AS ?out_col_name",
                                         summary = DBI::SQL(ifelse(zoi_metric == "cumulative", "sum", "max")),
                                         radius = DBI::SQL(radius[i]),
                                         zoi_limit = DBI::SQL(zoi_limit),
                                         input_geom = DBI::SQL(input_geom),
                                         infra_geom = DBI::SQL(infra_geom),
                                         out_col_name = DBI::SQL(output_column_name)[i])
          # dd_func <- "coalesce(exp(- log(1/?zoi_limit)/?radius * ST_Distance(pts.?input_geom, infra.?infra_geom)), 0)"
          # radius_div <- ifelse(output_type == "cumulative_zoi", 1, radius)
        }
      }
    }
    func_list[[i]] <- dd_func
  }

  funcs <- paste(func_list, collapse = ", ")

  base_query <- paste0("
SELECT pts.?input_points_id, ",
    funcs,
" FROM ?input_pts AS pts
    LEFT JOIN ?infra_layer AS infra
    ON ST_DWithin(pts.?input_geom, infra.?infra_geom, ?radius) ?another_condition
  GROUP BY pts.?input_points_id
limit ?lim;
")

  if(!is.null(output_table)) {
    create_tab <- DBI::sqlInterpolate(con, "CREATE TABLE IF NOT EXISTS ?out_name AS ", out_name = output_table)
    base_query <- paste0(create_tab, base_query)
  }

  qq <- DBI::sqlInterpolate(con, base_query,
                            # dist_decay_func = DBI::SQL(),
                            input_pts = DBI::SQL(input_points),
                            infra_layer = DBI::SQL(infrastructure_layer),
                            radius = DBI::SQL(radius),
                            input_geom = DBI::SQL(input_geom),
                            infra_geom = DBI::SQL(infra_geom),
                            input_points_id = DBI::SQL(input_id),
                            lim = DBI::SQL(limit),
                            another_condition = DBI::SQL(condition))
  qq
  if(verbose) print(qq)

  DBI::dbGetQuery(con, qq)
}
