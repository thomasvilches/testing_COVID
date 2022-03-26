using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(512), N=16, topology=:master_worker, exeflags = "--project=.")
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("covid19abm.jl")
@everywhere const cv=covid19abm


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
   
    # will return 6 dataframes. 1 total, 4 age-specific 
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
    end      

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors     

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    working = vcat([cdr[i].work for i = 1:nsims]...)
   
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g6 for i = 1:nsims]...)
    ag7 = vcat([cdr[i].g7 for i = 1:nsims]...) 

    mydfs = Dict("all" => allag, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6,"ag7" => ag7, "working"=>working)
    #mydfs = Dict("all" => allag, "working"=>working, "kids"=>kids)
    #mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    #c1 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :PRE, :MILD,:IISO, :HOS, :ICU, :DED), :_PREV)
    
    c1 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2,:LAT3, :HOS3, :ICU3, :DED3,:LAT4, :HOS4, :ICU4, :DED4,:LAT5, :HOS5, :ICU5, :DED5,:LAT6, :HOS6, :ICU6, :DED6), :_INC)
    c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2,:LAT3, :HOS3, :ICU3, :DED3,:LAT4, :HOS4, :ICU4, :DED4,:LAT5, :HOS5, :ICU5, :DED5,:LAT6, :HOS6, :ICU6, :DED6), :_PREV)
    
    #c2 = Symbol.((:LAT, :HOS, :ICU, :DED,:LAT2, :HOS2, :ICU2, :DED2), :_PREV)
    for (k, df) in mydfs
        println("saving dataframe sim level: $k")
        # simulation level, save file per health status, per age group
        #for c in vcat(c1..., c2...)
        for c in vcat(c1...)
        #for c in vcat(c2...)
            udf = unstack(df, :time, :sim, c) 
            fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
            CSV.write(fn, udf)
        end
        println("saving dataframe time level: $k")
        # time level, save file per age group
        #yaf = compute_yearly_average(df)       
        #fn = string("$(folderprefix)/timelevel_", k, ".dat")   
        #CSV.write(fn, yaf)       
    end
    

    writedlm(string(folderprefix,"/year_of_death.dat"),hcat([cdr[i].vector_dead for i=1:nsims]...))
    writedlm(string(folderprefix,"/npcr.dat"),hcat([cdr[i].npcr for i=1:nsims]...))
    writedlm(string(folderprefix,"/nra.dat"),hcat([cdr[i].nra for i=1:nsims]...))

    return mydfs
end


function create_folder(ip::cv.ModelParameters,province="ontario")
    
    #RF = string("heatmap/results_prob_","$(replace(string(ip.β), "." => "_"))","_vac_","$(replace(string(ip.vaccine_ef), "." => "_"))","_herd_immu_","$(ip.herd)","_$strategy","cov_$(replace(string(ip.cov_val)))") ## 
    main_folder = "/data/thomas-covid/testing_canada"
    #main_folder = "."
    
    RF = string(main_folder,"/results_prob_","$(replace(string(ip.β), "." => "_"))","_herd_immu_","$(ip.herd)","_idx_$(ip.file_index)_$(province)_strain_$(ip.strain)_scen_$(ip.scenariotest)_test_$(ip.test_ra)_eb_$(ip.extra_booster)_size_$(ip.size_threshold)") ##  
    
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end



function run_param_scen_cal(b::Float64,province::String="ontario",h_i::Int64 = 0,ic1::Int64=1,strains::Int64 = 1,index::Int64 = 0,scen::Int64 = 0,tra::Int64 = 0,eb::Int64 = 0,wpt::Int64 = 100,wpc::Float64 = 0.0,dayst::Vector{Int64} = [1;4],rc=[1.0],dc=[1],mt::Int64=300,vac::Bool=true,nsims::Int64=500)
    
    
    @everywhere ip = cv.ModelParameters(β=$b,fsevere = 1.0,fmild = 1.0,vaccinating = $vac,
    herd = $(h_i),start_several_inf=true,
    initialinf = $ic1,
    file_index = $index,
    modeltime=$mt, prov = Symbol($province),
    time_change_contact = $dc,
    change_rate_values = $rc,
    n_boosts = 1,
    scenariotest = $scen,
    extra_booster = $eb,
    size_threshold = $wpt,
    test_ra = $tra,
    testing_days = $dayst,
    strain = $strains,
    proportion_contacts_workplace = $wpc)

    folder = create_folder(ip,province)

    run(ip,nsims,folder)
   
end