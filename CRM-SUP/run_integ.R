lb <- c(qnorm(0.01, mu, sigma), qgamma(0.01, a_t, b_t), qgamma(0.01, a_g, b_g))
ub <- c(qnorm(0.99, mu, sigma), qgamma(0.99, a_t, b_t), qgamma(0.99, a_g, b_g))

the_denom <- integral3(denom, lb[1], ub[1], lb[2], ub[2], lb[3], ub[3], 
                       Yout=Yout, Dout=Dout, Wout=Wout, Zout=Zout, Nout=Nout,
                       mu=mu, sigma=sigma, a_t=a_t, b_t=b_t, a_g=a_g, b_g=b_g, reltol=0.001)

p1_post <- p2_post <- p3_post <- NULL
for (sched in 1:nsched)
{
  pp3 <- integral3(numer_p3, lb[1], ub[1], lb[2], ub[2], lb[3], ub[3], 
                  Yout=Yout, Dout=Dout, Wout=Wout, Zout=Zout, Nout=Nout,
                  mu=mu, sigma=sigma, a_t=a_t, b_t=b_t, a_g=a_g, b_g=b_g, 
                  doseval=skel[sched,3], reltol=0.001)/the_denom
  p3_post <- c(p3_post, pp3)
}
p_post <- cbind(p1_post, p2_post, p3_post)
