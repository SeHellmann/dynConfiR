fits <- structure(list(model = c("dynWEV", "PCRMt", "dynWEV", "PCRMt"),
                       sbj = c(1, 1, 8, 8),
                       negLogLik = c(20.6368391020115, 31.2510706889092, 17.5490953740148, 26.335957708293),
                       N = c(17, 17, 17, 17),
                       k = c(23, 17, 23, 19),
                       BIC = c(106.437585117316, 110.666768226774, 100.262097661323, 106.502968953654),
                       AICc = c(-57.2977503674057, -447.497858622182, -63.473237823399, -137.328084583414), AIC = c(87.2736782040229,
                                                                                                                    96.5021413778184, 81.0981907480295, 90.671915416586),
                       t0 = c(0, 0, 0.319525056128088, 0.311195873979081),
                       st0 = c(1.12963012317701, 1.38876864436006, 0.530550462493294, 0.550715158444022),
                       v1 = c(0.172998967009845, 0.803225303592474, 0.200143716989775, 0.103200570174599),
                       v2 = c(0.241457085865768, 0.257772690122886, 0.340016451531057, 0.173222271497872),
                       v3 = c(0.676480108664119,       0.612117874357785, 1.24517397378776, 0.665734848622986),
                       v4 = c(1.5499032887235,      1.72173493421176, 2.82777200274886, 1.6556441561065),
                       v5 = c(2.53678081878952,   3.39752560799776, 5.01116887971942, 4.009703275535),
                       thetaLower1 = c(1.17622428639492,    0.17179837861488, 0.789706736239607, 0.879092614245285),
                       thetaLower2 = c(1.17622428639492,   0.17179837861488, 1.42924978336366, 1.45498757702495),
                       thetaLower3 = c(1.17622428639492,      1.179837861488, 2.06880417585998, 2.02835133872529),
                       thetaLower4 = c(2.0508786601119,  1.82531013984407, 2.70909399201808, 2.58177942462658),
                       thetaUpper1 = c(1.20151430505538, 2.5445464439953, 0.789912481680675, 1.51912419362987),
                       thetaUpper2 = c(1.20151430505538, 3.11445464439953, 1.42982495424694, 2.07241846025026),
                       thetaUpper3 = c(1.20151430505538,      4.11445464439953, 2.06973742976804, 2.63297479420376),
                       thetaUpper4 = c(2.08163969688621,         5.06957951511996, 2.72942818864999, 3.21545615944865),
                       wrt = c(NA,        0.0362270956364099, NA, 0.0117852907511687),
                       wint = c(NA, 0.883174518432851,      NA, 0.872504068254756),
                       wx = c(NA, 0.109102515307999, NA, 0.124710987685301    ),
                       b = c(NA, 0.822557962816407, NA, 0.639762190335742),
                       a = c(1.76809729183024,    0.358796422821696, 1.33047749824592, 0.454392270307959),
                       z = c(0.144642348507182,        NA, 0.301234205760964, NA),
                       sz = c(0.117511684791096, NA, 0.120202480392568,   NA),
                       sigvis = c(0.228806568219654, NA, 0.189970554333006, NA),
                       tau = c(1.11663576818929,      NA, 0.989604987185369, NA),
                       w = c(0.296698197166943, NA, 0.489995452399694,    NA),
                       sig = c(0.231412989891031, NA, 0.49023842530622, NA),
                       svis = c(0.0321073472766225,   NA, 0.989920023402767, NA)),
                  row.names = c(NA, -4L), class = "data.frame")

preds <- predictConfModels(fits, subdivisions = 100, simult_conf = FALSE)
test_that("Prediction sums to 1", {
  expect_equal(sum(preds$p), 40, tolerance = 0.005)

})



# Regression coverage for optional parameters ---------------------------------

test_that("model-specific predictors tolerate missing optional parameters", {
  wev_params <- data.frame(
    a = 1, v = 0.3, t0 = 0.2, z = 0.5, tau = 0.1, w = 0.5,
    svis = 0.1, sigvis = 0.05
  )
  expect_silent(
    predictWEV_Conf(
      wev_params, model = "dynWEV", maxrt = 1, subdivisions = 10,
      simult_conf = FALSE, stop.on.error = FALSE, precision = 2,
      .progress = FALSE
    )
  )
  expect_silent(
    predictWEV_RT(
      wev_params, model = "dynWEV", maxrt = 1, subdivisions = 10,
      minrt = wev_params$t0, simult_conf = FALSE, scaled = FALSE,
      DistConf = NULL, precision = 2, .progress = FALSE
    )
  )

  rm_params <- data.frame(a = 1, b = 1.2, v = 0.4, t0 = 0.2, theta1 = 0.5)
  expect_silent(
    predictRM_Conf(
      rm_params, model = "IRM", maxrt = 1, subdivisions = 10,
      stop.on.error = FALSE, .progress = FALSE
    )
  )
  expect_silent(
    predictRM_RT(
      rm_params, model = "IRM", maxrt = 1, subdivisions = 10,
      minrt = rm_params$t0, scaled = FALSE, DistConf = NULL,
      .progress = FALSE
    )
  )

  mtl_params <- data.frame(
    v = 0.4, t0 = 0.2,
    mu_d1 = 0, mu_d2 = 0,
    s_v1 = 1, s_v2 = 1,
    s_d1 = 1, s_d2 = 1,
    rho_v = 0, rho_d = 0
  )
  expect_silent(
    predictMTLNR_Conf(
      mtl_params, maxrt = 1, subdivisions = 10,
      stop.on.error = FALSE, .progress = FALSE
    )
  )
  expect_silent(
    predictMTLNR_RT(
      mtl_params, maxrt = 1, subdivisions = 10,
      minrt = mtl_params$t0, scaled = FALSE, DistConf = NULL,
      .progress = FALSE
    )
  )

  dd_params <- data.frame(a = 1, v = 0.2, t0 = 0.2, z = 0.5, sz = NA_real_, sv = NA_real_)
  expect_silent(
    predictDDConf_Conf(
      dd_params, maxrt = 1, subdivisions = 10,
      stop.on.error = FALSE, .progress = FALSE
    )
  )
  expect_silent(
    predictDDConf_RT(
      dd_params, maxrt = 1, subdivisions = 10,
      minrt = dd_params$t0, scaled = FALSE, DistConf = NULL,
      .progress = FALSE
    )
  )
})


test_that("predictConfModels and predictRTModels keep identifiers leading", {
  params <- data.frame(
    model = "dynWEV", sbj = 1,
    a = 1, v = 0.4, t0 = 0.2, z = 0.5, tau = 0.1,
    w = 0.5, svis = 0.1, sigvis = 0.05, theta1 = 0.6
  )
  conf_dist <- predictConfModels(
    params, maxrt = 1, subdivisions = 5,
    simult_conf = FALSE, stop.on.error = FALSE,
    .progress = FALSE, parallel = FALSE
  )
  expect_identical(names(conf_dist)[1:2], c("sbj", "model"))

  rt_dist <- predictRTModels(
    params, maxrt = 1, subdivisions = 5, minrt = params$t0,
    simult_conf = FALSE, scaled = FALSE, DistConf = NULL,
    .progress = FALSE, parallel = FALSE
  )
  expect_identical(names(rt_dist)[1:2], c("sbj", "model"))
})


#as.matrix(jobs)[2,]
