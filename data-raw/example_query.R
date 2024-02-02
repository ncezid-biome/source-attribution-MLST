devtools::load_all()

li_path <- system.file("extdata/example_query.csv", package = "sourcerer")

example_query <- readr::read_csv(li_path)

usethis::use_data(example_query, overwrite = TRUE)
