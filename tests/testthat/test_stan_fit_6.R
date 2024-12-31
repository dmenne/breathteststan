test_that("Bad records are skipped when multiple in one file ", {
  skip_on_cran()
  library(breathtestcore)
  chains = 1
  student_t_df = 10
  dose = 100
  iter = 100
  # This file has two bad records
  filename = btcore_file("short_record.xml")
  xml_data = read_any_breathtest(filename)
  data = cleanup_data(xml_data)
  fit = stan_fit(data, dose = dose, student_t_df = student_t_df,
                 chains = chains, iter = iter  )
  expect_s3_class(fit, "breathtestfit")
})

