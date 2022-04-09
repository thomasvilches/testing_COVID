#for the paper

#run_param_scen_cal(b,province,h_i,ic1,strains,index,scen,tra,eb,wpt,dayst,rc,dc,mt,vac,nsims)
#run_param_scen_cal(0.042,"ontario",10,1,3,2,0,0,0,50,[1],1.0,0.0,50)
#run_param_scen_cal(0.042,"ontario",10,1,3,1,0,0,0,50,[1],1.0,0.0,365)
#0.35 for omicron generates a R0 of 0.84
#= 
#run_param_scen_cal(b,province,h_i,ic1,strains,index,scen,tra,eb,wpt,dayst,rc,dc,mt,vac,nsims)
run_param_scen_cal(0.042,"ontario",10,1,3,1,0,0,0,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,1,0,0,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,2,2,0,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,3,2,0,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,4,2,0,50,[1,4],1.0,0.0,365,1,365)
##
run_param_scen_cal(0.042,"ontario",10,1,3,2,0,0,0,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,1,0,0,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,2,2,0,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,3,2,0,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,4,2,0,50,[1,4],1.0,0.0,365,1,112)

run_param_scen_cal(0.042,"ontario",10,1,3,3,0,0,0,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,1,0,0,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,2,2,0,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,3,2,0,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,4,2,0,50,[1,4],1.0,0.0,365,1,224)


run_param_scen_cal(0.042,"ontario",10,1,3,4,0,0,0,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,1,0,0,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,2,2,0,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,3,2,0,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,4,2,0,50,[1,4],1.0,0.0,365,112,112)


run_param_scen_cal(0.042,"ontario",10,1,3,5,0,0,0,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,1,0,0,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,2,2,0,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,3,2,0,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,4,2,0,50,[1,4],1.0,0.0,365,112,253)

 =#

#################################################
########### 20% booster
#################################################

#run_param_scen_cal(b,province,h_i,ic1,strains,index,scen,tra,eb,wpt,dayst,rc,dc,mt,vac,nsims)
run_param_scen_cal(0.042,"ontario",10,1,3,1,0,0,20,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,1,0,20,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,2,2,20,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,3,2,20,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,4,2,20,50,[1,4],1.0,0.0,365,1,365)
##
run_param_scen_cal(0.042,"ontario",10,1,3,2,0,0,20,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,1,0,20,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,2,2,20,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,3,2,20,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,4,2,20,50,[1,4],1.0,0.0,365,1,112)

run_param_scen_cal(0.042,"ontario",10,1,3,3,0,0,20,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,1,0,20,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,2,2,20,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,3,2,20,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,4,2,20,50,[1,4],1.0,0.0,365,1,224)


run_param_scen_cal(0.042,"ontario",10,1,3,4,0,0,20,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,1,0,20,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,2,2,20,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,3,2,20,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,4,2,20,50,[1,4],1.0,0.0,365,112,112)


run_param_scen_cal(0.042,"ontario",10,1,3,5,0,0,20,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,1,0,20,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,2,2,20,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,3,2,20,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,4,2,20,50,[1,4],1.0,0.0,365,112,253)




#################################################
########### 80% booster
#################################################

#run_param_scen_cal(b,province,h_i,ic1,strains,index,scen,tra,eb,wpt,dayst,rc,dc,mt,vac,nsims)
run_param_scen_cal(0.042,"ontario",10,1,3,1,0,0,80,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,1,0,80,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,2,2,80,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,3,2,80,50,[1,4],1.0,0.0,365,1,365)
run_param_scen_cal(0.042,"ontario",10,1,3,1,4,2,80,50,[1,4],1.0,0.0,365,1,365)
##
run_param_scen_cal(0.042,"ontario",10,1,3,2,0,0,80,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,1,0,80,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,2,2,80,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,3,2,80,50,[1,4],1.0,0.0,365,1,112)
run_param_scen_cal(0.042,"ontario",10,1,3,2,4,2,80,50,[1,4],1.0,0.0,365,1,112)

run_param_scen_cal(0.042,"ontario",10,1,3,3,0,0,80,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,1,0,80,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,2,2,80,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,3,2,80,50,[1,4],1.0,0.0,365,1,224)
run_param_scen_cal(0.042,"ontario",10,1,3,3,4,2,80,50,[1,4],1.0,0.0,365,1,224)


run_param_scen_cal(0.042,"ontario",10,1,3,4,0,0,80,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,1,0,80,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,2,2,80,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,3,2,80,50,[1,4],1.0,0.0,365,112,112)
run_param_scen_cal(0.042,"ontario",10,1,3,4,4,2,80,50,[1,4],1.0,0.0,365,112,112)


run_param_scen_cal(0.042,"ontario",10,1,3,5,0,0,80,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,1,0,80,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,2,2,80,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,3,2,80,50,[1,4],1.0,0.0,365,112,253)
run_param_scen_cal(0.042,"ontario",10,1,3,5,4,2,80,50,[1,4],1.0,0.0,365,112,253)






