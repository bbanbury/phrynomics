test_that("Calculate Missing Data"){
#Test for making sure that data is being read in using ReadSNP and that missing data is being calculated correctly. 
  data(fakeData)
  snpfakeData <- ReadSNP(fakeData)
  expect_identical(CalculateMissingData(fakeData), CalculateMissingData(snpfakeData))
}