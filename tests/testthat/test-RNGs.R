test_that("r2DSD works", {
  expect_true(all(r2DSD(200, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2)$response %in% c(-1,0,1)))
  expect_true(all(r2DSD(200, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2)$rt >= 0))

  expect_length(r2DSD(20, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2)$rt, 20)
  expect_true(any(r2DSD(20, tau=1, a=30, v=0.1, t0=0, z =0.5, sv=0.2, maxrt = 2)$response==0))
  expect_error(r2DSD(20, tau=1, a=2, v=0.1, t0=-1, z =0.5, sv=0.2, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(r2DSD(20, tau=1, a=2, v=0.1, t0=-1, z =0.5, sv=0.2, maxrt = 5, stop_on_error = FALSE)== 0))

  expect_error(r2DSD(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=-0.2, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(r2DSD(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=-0.2, maxrt = 5, stop_on_error = FALSE)== 0))

  expect_error(r2DSD(20, tau=1, a=2, v=0.1, t0=1, z =-0.5, sv=0.2, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(r2DSD(20, tau=1, a=2, v=0.1, t0=1, z =1.5, sv=0.2, maxrt = 5, stop_on_error = FALSE)== 0))
})

test_that("rdynaViTE works", {
  expect_true(all(rdynaViTE(200, tau=1, a=2, v=c(0.5,-0.5), t0=0, z =0.5, sv=0.2, w= 0.8)$response %in% c(-1,0,1)))
  expect_true(all(rdynaViTE(200, tau=1, a=2, v=0.5, t0=0, z =c(0.4,0.5, 0.6), sv=0.2, w= 0.8)$rt >= 0))

  expect_length(rdynaViTE(20, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, w= 0.8)$rt, 20)
  expect_true(any(rdynaViTE(20, tau=1, a=30, v=0.1, t0=0, z =0.5, sv=0.2, maxrt = 2, w= 0.8)$response==0))
  expect_error(rdynaViTE(20, tau=1, a=2, v=0.1, t0=-1, z =0.5, sv=0.2, w= 0.8, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(rdynaViTE(20, tau=1, a=2, v=0.1, t0=-1, z =0.5, sv=0.2, w= 0.8, maxrt = 5, stop_on_error = FALSE)== 0))

  expect_error(rdynaViTE(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=0.2, svis =0, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(rdynaViTE(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=0.2, svis = -1, maxrt = 5, stop_on_error = FALSE)== 0))

  expect_error(rdynaViTE(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=0.2, w= 1.2, maxrt = 5, stop_on_error = TRUE))
  expect_true(all(rdynaViTE(20, tau=1, a=2, v=0.1, t0=1, z =0.5, sv=0.2, w= -0.8, maxrt = 5, stop_on_error = FALSE)== 0))
})



test_that("rPCRM works", {
  expect_true(all(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$response %in% c(0,1,2)))
  expect_true(all(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$rt >= 0))

  expect_length(rPCRM(20, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                      wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$rt, 20)
  expect_true(any(rPCRM(250, mu1=1.5, mu2=1.8, a=25, b=25,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1, maxrt = 15)$response==0))
  expect_error(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=-1, st0=0.2, s=1))
  expect_error(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=-0.2, s=1))
  expect_error(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=-1))
  expect_error(rPCRM(250, mu1=1.5, mu2=1.8, a=-2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=1))
  expect_error(rPCRM(250, mu1=1.5, mu2=1.8, a=2.5, b=-2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=1))
})


test_that("rIRM works", {
  expect_true(all(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$response %in% c(0,1,2)))
  expect_true(all(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$rt >= 0))

  expect_length(rIRM(20, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                      wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1)$rt, 20)
  expect_true(any(rIRM(250, mu1=1.5, mu2=1.8, a=25, b=25,
                        wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1, maxrt = 15)$response==0))
  expect_error(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=-1, st0=0.2, s=1))
  expect_error(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=-0.2, s=1))
  expect_error(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=-1))
  expect_error(rIRM(250, mu1=1.5, mu2=1.8, a=-2.5, b=2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=1))
  expect_error(rIRM(250, mu1=1.5, mu2=1.8, a=2.5, b=-2.5,
                     wx=0.5, wrt=0.1, wint=0.4, t0=1, st0=0.2, s=1))
})
