test_that("2DSD works", {
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", 0.4, 2.5, a=2, v=0.5, t0=0, z =0.5, sv=0.2, tau=1),
               c(0, 0.127923586525416, 0.0668243025079862, 0.0342270990259822,
                 0.0175585010435632, 0.0090222312746275, 0.00464316676344542))
  ## test for z (relative) out of range
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", 0.4, 2.5, a=2, v=0.5, t0=0, z =1.1, sv=0.2, tau=1, s=1, stop_on_error = FALSE),
               rep(0, 7))
  ## test for absolute z working (if as relative it would be out of range)
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", 0.4, 2.5, a=2, v=0.5, t0=0, z =1.1, sv=0.2, tau=1, s=1, stop_on_error = TRUE, z_absolute = TRUE),
               c(0, 0.115664728481645, 0.0659737402126581, 0.0342454229147437,
                 0.0176076118813728, 0.00905219904467997, 0.00465973138080759))
  ## test for relative z giving same results as transformed absolute z
  expect_equal(d2DSD(seq(0, 3, by=0.5),  "lower", 0.4, 2.5, a=2, v=0.5, t0=0, z =0.5, sv=0.2,tau=1),
               d2DSD(seq(0, 3, by=0.5), "lower", 0.4, 2.5, a=2, v=0.5, t0=0, z =1, sv=0.2, tau=1, s=1, z_absolute = TRUE))
  ## test for effect of t0
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", -100, 2.5, a=2, v=0.5, t0=0, z =0.2, sv=0.2, tau=1, s=1, stop_on_error = TRUE),
               d2DSD(seq(0, 3, by=0.5)+0.3, "lower", -100, 2.5, a=2, v=0.5, t0=0.3, z =0.2, sv=0.2, tau=1, s=1, stop_on_error = TRUE))
  ## test for scaling effect of s (diffusion constant)
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", 0.5*2, 2.5*2, a=2*2, v=0.5*2, t0=0, z =0.2, sv=0.2*2, tau=1, s=1*2, stop_on_error = TRUE),
               d2DSD(seq(0, 3, by=0.5), "lower", 0.5, 2.5, a=2, v=0.5, t0=0, z =0.2, sv=0.2, tau=1, s=1, stop_on_error = TRUE))
  ## test for response symmetry
  expect_equal(d2DSD(seq(0, 3, by=0.5), "lower", -1, 1, a=2, v=0.5, t0=0, z =0.2, sv=0.2, tau=1, s=1, stop_on_error = TRUE),
               d2DSD(seq(0, 3, by=0.5), "upper", -1, 1,  a=2, v=-0.5, t0=0, z =1-0.2, sv=0.2, tau=1, s=1, stop_on_error = TRUE))
  ## test for equivalence wrt to lambda
  expect_equal(d2DSD(c(0.1, .3,.7,1.2), "upper", 0.5 /((c(0.1, .3,.7,1.2)+2)), 1/((c(0.1, .3,.7,1.2)+2)),  2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2,1),
               d2DSD(c(0.1, .3,.7,1.2), "upper", 0.5 /sqrt((c(0.1, .3,.7,1.2)+2)), 1/sqrt((c(0.1, .3,.7,1.2)+2)),  2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2,0.5))
  expect_equal(d2DSD(c(0.1, .3,.7,1.2), "lower", 0.5, 1, 2,  0.5,  0,0.5, 0, 0, 0.3, 0, 2), # default lambda = 0
               d2DSD(c(0.1, .3,.7,1.2), "lower", 0.5 /(c(0.1, .3,.7,1.2)+2)^2, 1/(c(0.1, .3,.7,1.2)+2)^2, 2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2,2))
})


test_that("dynWEV works", {
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, w=0.5, sigvis=0.2, svis=1, s=1),
               c(0, 0.147064426866004, 0.0798344162131249, 0.0414020207786117,
                 0.0210877315338026, 0.0106006919740067, 0.00527825242394075))
  ## test for scaling effect of s (diffusion constant)
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", 0.4*2, 2.5*2, tau=1, a=2*2, v=0.5*2, t0=0, z =0.5, sv=0.2*2, w=0.5, sigvis=0.2*2, svis=1*2, s=1*2, precision = 10),
               dWEV(seq(0, 3, by=0.5), "lower", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, w=0.5, sigvis=0.2, svis=1, s=1, precision = 10))
  ## test for effect of t0
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, w=0.5, sigvis=0.2, svis=1, s=1),
               dWEV(seq(0, 3, by=0.5)+0.4, "lower", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0.4, z =0.5, sv=0.2, w=0.5, sigvis=0.2, svis=1, s=1))
  ## test with negative lower confidence threshold and different parameters
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", -100, 100, tau=1, a=2, v=-0.5, t0=0, z =0.2, sv=0.2, w=0.001, sigvis=0.4, svis=1, s=1),
               c(0, 0.43714364613092, 0.153514811046022, 0.0724144673283747,
                 0.0361455558387259, 0.0182440173062963, 0.00924134169763999))
  ## test with explicitly giving parameter muvis
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", -100, 100, tau=1, a=2, v=-0.5, t0=0, z =0.2, sv=0.2, w=0.001, sigvis=0.4, svis=1, s=1),
               dWEV(seq(0, 3, by=0.5), "lower", -100, 100, tau=1, a=2, v=-0.5, t0=0, z =0.2, sv=0.2, w=0.001, muvis = 0.5, sigvis=0.4, svis=1, s=1))
  ## test for equivalence wrt to lambda
  expect_equal(dWEV(c(0.1, .3,.7,1.2), "upper", 0.5 /((c(0.1, .3,.7,1.2)+2)), 1/((c(0.1, .3,.7,1.2)+2)), 2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2, 0.5, lambda=1),
               dWEV(c(0.1, .3,.7,1.2), "upper", 0.5 /sqrt((c(0.1, .3,.7,1.2)+2)), 1/sqrt((c(0.1, .3,.7,1.2)+2)),  2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2, 0.5, lambda=0.5))
  expect_equal(dWEV(c(0.1, .3,.7,1.2), "lower", 0.5, 1, 2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2, 0.5), # default lambda = 0
               dWEV(c(0.1, .3,.7,1.2), "lower", 0.5 /(c(0.1, .3,.7,1.2)+2)^2, 1/(c(0.1, .3,.7,1.2)+2)^2, 2, 0.5,  0,0.5, 0, 0, 0.3, 0, 2,  0.5, lambda=2))
})

test_that("dynWEV with w=0 equals 2DSD", {
  expect_equal(dWEV(seq(0, 3, by=0.5), "upper", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, st0=0, sz=0,w=1-.Machine$double.eps, sigvis=0.2, svis=1, s=1, precision = 3),
               d2DSD(seq(0, 3, by=0.5), "upper", 0.4, 2.5, tau=1, a=2, v=0.5, t0=0, z =0.5, sv=0.2, s=1, st0=0, sz=0, precision = 3))
  expect_equal(dWEV(seq(0, 3, by=0.5), "lower", 0.4, 2.5, tau=1, a=2, v=-0.5, t0=0, z =0.5, sv=0.2, st0=0, sz=0,w=1-.Machine$double.eps, sigvis=0.2, svis=1, s=1, precision = 3),
               d2DSD(seq(0, 3, by=0.5), "lower", 0.4, 2.5, tau=1, a=2, v=-0.5, t0=0, z =0.5, sv=0.2, s=1, st0=0, sz=0, precision = 3))
})





test_that("IRM works", {
  expect_equal(dIRM(seq(0, 3, by=0.5), response=1, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1),
               c(0, 0.0476788999581483, 0.25567473708014, 0.0990072008828482,
                 0.0236905178296517, 0.00498839744174124, 0.00101269796804896))
  ## Test for scaling of s
  expect_equal(dIRM(seq(0, 3, by=0.5), response=1, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1),
               dIRM(seq(0, 3, by=0.5), response=1, mu1=1.5*2, mu2=1.8*2, a=2.5*2, b=2.5*2,
                    th1=1*2, th2=10*2, wx=0.5, wrt=0.1*2, wint=0.4, t0=0, st0=0.2, s=1*2))
  ## Test for RT-shift with t0 using different parameters
  expect_equal(dIRM(seq(0, 3, by=0.5)+0.4, response=2, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8),
               dIRM(seq(0, 3, by=0.5), response=2, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1.8))
  ## Test for failing, if parameter set is out of range
  expect_error(dIRM(1.5, response=1, mu1=0.5, mu2=0.8, a=-2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8), "thresholds")
  expect_error(dIRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=-2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8), "thresholds")
  expect_error(dIRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=-0.4, st0=0.2, s=1.8), "t0")
  expect_error(dIRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=-0.2, s=1.8), "st0")
  expect_error(dIRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=-1.8), "positive")
})


test_that("PCRM works", {
  expect_equal(dPCRM(seq(0, 3, by=0.5), response=1, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1),
               c(0, 0.0519309447177127, 0.396833386300028, 0.107050177325627,
                 0.0105874394413118, 0.000729360533047754, 4.3009284047349e-05))
  ## Test for scaling of s
  expect_equal(dPCRM(seq(0, 3, by=0.5), response=1, mu1=1.5, mu2=1.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1),
               dPCRM(seq(0, 3, by=0.5), response=1, mu1=1.5*2, mu2=1.8*2, a=2.5*2, b=2.5*2,
                    th1=1*2, th2=10*2, wx=0.5, wrt=0.1*2, wint=0.4, t0=0, st0=0.2, s=1*2))
  ## Test for RT-shift with t0 using different parameters
  expect_equal(dPCRM(seq(0, 3, by=0.5)+0.4, response=2, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8),
               dPCRM(seq(0, 3, by=0.5), response=2, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0, st0=0.2, s=1.8))
  ## Test for failing, if parameter set is out of range
  expect_error(dPCRM(1.5, response=1, mu1=0.5, mu2=0.8, a=-2.5, b=2.5,
                    th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8), "thresholds")
  expect_error(dPCRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=-2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=1.8), "thresholds")
  expect_error(dPCRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=-0.4, st0=0.2, s=1.8), "t0")
  expect_error(dPCRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=-0.2, s=1.8), "st0")
  expect_error(dPCRM(1.5, response=1, mu1=0.5, mu2=0.8, a=2.5, b=2.5,
                     th1=1, th2=10, wx=0.5, wrt=0.1, wint=0.4, t0=0.4, st0=0.2, s=-1.8), "positive")
})

