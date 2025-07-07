autoplot(fs_comp_trans[[1]], "PE.A") +
  ggtitle("APC-A Histogram") +
  xlab("APC-A")+
  ylab("Count")

autoplot(fs_noncomp_trans[[7]], "APC.A") +
  ggtitle("APC-A Histogram") +
  xlab("APC-A")+
  ylab("Count")
  
autoplot(fs_controls_noncomp_trans[[3]], "PE.A") +
  ggtitle("APC-A Histogram") +
  xlab("APC-A")+
  ylab("Count")

fs_controls_noncomp_trans[[3]]
