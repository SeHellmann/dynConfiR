library(dynConfiR)

curve(d2DSD(x, "upper", -Inf, Inf, tau=1, a=2, v=0.4, sz=0.2, sv=0.9), xlim=c(0, 2), lty=2)
curve(d2DSD(x, "lower", -Inf, Inf, tau=1, a=2, v=0.4, sz=0.2, sv=0.9), col="red", lty=2, add=TRUE)
curve(d2DSD(x, "upper", -Inf, Inf, tau=1, a=2, v=0.4),add=TRUE)
curve(d2DSD(x, "lower", -Inf, Inf, tau=1, a=2, v=0.4), col="red", add=TRUE)
# Generate a random sample
dfu <- r2DSD(5000, a=2,v=0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1)
# Same RT distribution but upper and lower responses changed
dfl <- r2DSD(50, a=2,v=-0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1)
head(dfu)
d2DSD(dfu, th1=-Inf, th2=Inf, a=2, v=.5)[1:5]
# Scaling diffusion parameters leads do same density values
s <- 2
d2DSD(dfu, th1=-Inf, th2=Inf, a=2*s, v=.5*s, s=2)[1:5]
if (requireNamespace("ggplot2", quietly = TRUE)) {
 require(ggplot2)
 ggplot(dfu, aes(x=rt, y=conf))+
   stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
   facet_wrap(~response)
}
boxplot(conf~response, data=dfu)
# Restricting to specific confidence region
dfu <- dfu[dfu$conf >0 & dfu$conf <1,]
d2DSD(dfu, th1=0, th2=1, a=2, v=0.5)[1:5]
# If lower confidence threshold is higher than the upper, the function throws an error,
# except when stop_on_error is FALSE
d2DSD(dfu[1:5,], th1=1, th2=0, a=2, v=0.5, stop_on_error = FALSE)








#
# user@kap-u18-pc3:~/dynConfiR$ R -d "valgrind --leak-check=full --track-origins=yes -s" -f run_examples.R
# ==79109== Memcheck, a memory error detector
# ==79109== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
# ==79109== Using Valgrind-3.18.1 and LibVEX; rerun with -h for copyright info
# ==79109== Command: /usr/lib/R/bin/exec/R -f run_examples.R
# ==79109==
#
# R version 4.1.2 (2021-11-01) -- "Bird Hippie"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)
#
# R ist freie Software und kommt OHNE JEGLICHE GARANTIE.
# Sie sind eingeladen, es unter bestimmten Bedingungen weiter zu verbreiten.
# Tippen Sie 'license()' or 'licence()' für Details dazu.
#
# R ist ein Gemeinschaftsprojekt mit vielen Beitragenden.
# Tippen Sie 'contributors()' für mehr Information und 'citation()',
# um zu erfahren, wie R oder R packages in Publikationen zitiert werden können.
#
# Tippen Sie 'demo()' für einige Demos, 'help()' für on-line Hilfe, oder
# 'help.start()' für eine HTML Browserschnittstelle zur Hilfe.
# Tippen Sie 'q()', um R zu verlassen.
#
# > library(dynConfiR)
# >
# > curve(d2DSD(x, "upper", -Inf, Inf, tau=1, a=2, v=0.4, sz=0.2, sv=0.9), xlim=c(0, 2), lty=2)
# > curve(d2DSD(x, "lower", -Inf, Inf, tau=1, a=2, v=0.4, sz=0.2, sv=0.9), col="red", lty=2, add=TRUE)
# > curve(d2DSD(x, "upper", -Inf, Inf, tau=1, a=2, v=0.4),add=TRUE)
# > curve(d2DSD(x, "lower", -Inf, Inf, tau=1, a=2, v=0.4), col="red", add=TRUE)
# > # Generate a random sample
# > dfu <- r2DSD(5000, a=2,v=0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1)
# > # Same RT distribution but upper and lower responses changed
# > dfl <- r2DSD(50, a=2,v=-0.5,t0=0,z=0.5,d=0,sz=0,sv=0, st0=0,  tau=1, s=1)
# > head(dfu)
#     rt response      conf
# 1 0.95        1 0.8846157
# 2 0.74        1 1.5253340
# 3 1.99        1 1.0620982
# 4 0.57        1 2.9496727
# 5 0.64        1 1.4961185
# 6 2.99        1 2.1624011
# > d2DSD(dfu, th1=-Inf, th2=Inf, a=2, v=.5)[1:5]
# [1] 0.35608557 0.47282760 0.08669331 0.59043553 0.53979489
# > # Scaling diffusion parameters leads do same density values
# > s <- 2
# > d2DSD(dfu, th1=-Inf, th2=Inf, a=2*s, v=.5*s, s=2)[1:5]
# [1] 0.35608557 0.47282760 0.08669331 0.59043553 0.53979489
# > if (requireNamespace("ggplot2", quietly = TRUE)) {
# +  require(ggplot2)
# +  ggplot(dfu, aes(x=rt, y=conf))+
# +    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
# +    facet_wrap(~response)
# + }
# Lade nötiges Paket: ggplot2
# > boxplot(conf~response, data=dfu)
# > # Restricting to specific confidence region
# > dfu <- dfu[dfu$conf >0 & dfu$conf <1,]
# > d2DSD(dfu, th1=0, th2=1, a=2, v=0.5)[1:5]
# [1] 0.08607669 0.06833802 0.08043123 0.07117993 0.16992300
# > # If lower confidence threshold is higher than the upper, the function throws an error,
# > # except when stop_on_error is FALSE
# > d2DSD(dfu[1:5,], th1=1, th2=0, a=2, v=0.5, stop_on_error = FALSE)
# [1] 0 0 0 0 0
# >
# ==79109==
# ==79109== HEAP SUMMARY:
# ==79109==     in use at exit: 110,805,665 bytes in 23,477 blocks
# ==79109==   total heap usage: 105,392 allocs, 81,915 frees, 331,030,247 bytes allocated
# ==79109==
# ==79109== LEAK SUMMARY:
# ==79109==    definitely lost: 0 bytes in 0 blocks
# ==79109==    indirectly lost: 0 bytes in 0 blocks
# ==79109==      possibly lost: 0 bytes in 0 blocks
# ==79109==    still reachable: 110,805,665 bytes in 23,477 blocks
# ==79109==                       of which reachable via heuristic:
# ==79109==                         newarray           : 4,264 bytes in 1 blocks
# ==79109==         suppressed: 0 bytes in 0 blocks
# ==79109== Reachable blocks (those to which a pointer was found) are not shown.
# ==79109== To see them, rerun with: --leak-check=full --show-leak-kinds=all
# ==79109==
# ==79109== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
#
