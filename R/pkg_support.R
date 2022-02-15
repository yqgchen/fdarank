# This file is a cheat to minimize the false positives flagged during R CMD check such as
# "no visible binding for global variable 'j'"
# when in smoRank(),
# ise <- foreach(j = idxInnerPts, .combine=c ) %doforeach% {
#   est <- cdfsm(tin=tin[-j], yin=yin[-j], tout=tin[j], yout=yWorkGrid,
#                ybw=ybwCand[ybwIdx], tbw=tbwCand[tbwIdx], H=H, K=K)
#   err <- ((yWorkGrid >= yin[j]) - est)^2
#   pracma::trapz(x = yWorkGrid, y = err)
# }
utils::globalVariables(c("j"))
