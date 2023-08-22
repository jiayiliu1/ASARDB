context("SelAgg")

load('~/Documents/Projects/R/asarDB/china.RData')

taxS<-'genus'
taxVal<-'Geobacter'
taxVaL<-c('Geobacter','Burkholderia')
funS<-'FUN3'
funVal<-'Isoprenoids'
samples<-'mgm4714679.3'
taxAgg<-'species'
funAgg<-'FUN2'

test_that('selagg perform well on atomic full input',{
  funtax<-selagg(funtaxall,
                 taxS=taxS,
                 taxVal=taxVal,
                 funS=funS,
                 funVal=funVal,
                 taxAgg=taxAgg,
                 funAgg=funAgg,
                 samples=samples)
  expect_is(funtax,'data.table')
  expect_equal(dim(funtax),c(14,3))
})


test_that('selagg perform well on full input with multiple samples',{
  samples<-c("mgm4714677.3","mgm4714679.3")
  funtax<-selagg(funtaxall,
                 taxS=taxS,
                 taxVal=taxVal,
                 funS=funS,
                 funVal=funVal,
                 taxAgg=taxAgg,
                 funAgg=funAgg,
                 samples=samples)
  expect_is(funtax,'data.table')
  expect_equal(dim(funtax),c(14,4))
})

test_that('selagg perform well on full input with multiple taxons',{
  funtax<-selagg(funtaxall,
                 taxS=taxS,
                 taxVal=taxVaL,
                 funS=funS,
                 funVal=funVal,
                 taxAgg=taxAgg,
                 funAgg=funAgg,
                 samples=samples)
  expect_is(funtax,'data.table')
  expect_equal(dim(funtax),c(14,4))
})

