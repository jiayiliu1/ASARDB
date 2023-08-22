context("Dimensions")

tax<-c('usp','species','genus','family','order','class','phylum','domain','kingdom')

test_that('dimNames returns all dimensions for taxonomy',{
  expect_type(dimNames('taxonomy'),'character')
  expect_true(length(dimNames('taxonomy'))>0)
  expect_false(any(is.na(match(tolower(dimNames('taxonomy')),tax))))
})

test_that('dimNames returns all dimensions for taxonomy with abbreviated call',{
  expect_false(any(is.na(match(tolower(dimNames('taxonomy')),tax))))
  expect_false(any(is.na(match(tolower(dimNames('tax')),tax))))
  expect_false(any(is.na(match(tolower(dimNames('t')),tax))))

})

fun<-c('ufun','fun1','fun2','fun3','fun4','ko')
test_that('dimNames returns all dimensions for functions',{
  expect_type(dimNames('functions'),'character')
  expect_true(length(dimNames('functions'))>0)
  expect_false(any(is.na(match(tolower(dimNames('functions')),fun))))
})

test_that('dimNames returns all dimensions for functions with abbreviated call',{
  expect_false(any(is.na(match(tolower(dimNames('functions')),fun))))
  expect_false(any(is.na(match(tolower(dimNames('fun')),fun))))
  expect_false(any(is.na(match(tolower(dimNames('f')),fun))))
})

test_that('dimNames throws an error if name neither function nor taxonomy',{
  expect_error(dimNames('samples'),regexp = "'arg' should be one of")
})

good_matrix<-data.frame(usp=c('test1', 'test2', 'test3'),
                        species=c('test1', 'test2', 'test3'),
                        genus=c('test1', 'test2', 'test3'),
                        family=c('test1', 'test2', 'test3'),
                        order=c('test1', 'test2', 'test3'),
                        class=c('test1', 'test2', 'test3'),
                        phylum=c('test1', 'test2', 'test3'),
                        domain=c('test1', 'test2', 'test3'),
                        md5=c('test1', 'test2', 'test3'),
                        ufun=c('test1', 'test2', 'test3'),
                        FUN2=c('test1', 'test2', 'test3'),
                        FUN3=c('test1', 'test2', 'test3'),
                        FUN4=c('test1', 'test2', 'test3'),
                        mgm4714659.3=c('test1', 'test2', 'test3'),
                        mgm4714661.3=c('test1', 'test2', 'test3'),
                        mgm4714663.3=c('test1', 'test2', 'test3'),
                        mgm4714665.3=c('test1', 'test2', 'test3'),
                        mgm4714667.3=c('test1', 'test2', 'test3'),
                        mgm4714669.3=c('test1', 'test2', 'test3'),
                        mgm4714671.3=c('test1', 'test2', 'test3'),
                        mgm4714673.3=c('test1', 'test2', 'test3'),
                        mgm4714675.3=c('test1', 'test2', 'test3'),
                        mgm4714677.3=c('test1', 'test2', 'test3'))
sn<-c("mgm4714659.3", "mgm4714661.3", "mgm4714663.3", "mgm4714665.3", "mgm4714667.3", "mgm4714669.3", "mgm4714671.3", "mgm4714673.3", "mgm4714675.3", "mgm4714677.3")

test_that("getSampleNames return all sample names",{
  expect_equivalent(getSampleNames(good_matrix),sn)
  ridx<-sample.int(dim(good_matrix)[2])
  expect_equivalent(sort(getSampleNames(good_matrix[,ridx])),sort(sn))

})

test_that('getDim return all dimensions for taxonomy',{
  idx<-getDim(good_matrix,'taxonomy')
  expect_equal(length(idx),8)
  expect_equivalent(idx,c(1:8))
  ridx<-sample.int(dim(good_matrix)[2])
  expect_equivalent(getDim(good_matrix[,ridx],'taxonomy'),match(1:8,ridx))
})
test_that('getDim return all dimensions for function',{
  idx<-getDim(good_matrix,'functions')
  expect_equal(length(idx),4)
  expect_equivalent(idx,c(10:13))
  ridx<-sample.int(dim(good_matrix)[2])
  expect_equivalent(getDim(good_matrix[,ridx],'functions'),match(10:13,ridx))
})
test_that('getDim return all dimensions for samples',{
  idx<-getDim(good_matrix,'samples')
  expect_equal(length(idx),10)
  expect_equivalent(idx,c(14:23))
  ridx<-sample.int(dim(good_matrix)[2])
  expect_equivalent(getDim(good_matrix[,ridx],'samples'),match(14:23,ridx))
})
# test_that('dimNames returns all dimensions for KEGG orthology',{
#
# })
