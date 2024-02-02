devtools::load_all()

li_path <- system.file("extdata/Listeria_isolates.csv.gz", package = "sourcerer")

Listeria_isolates <- readr::read_csv(li_path)

usethis::use_data(Listeria_isolates, overwrite = TRUE)
