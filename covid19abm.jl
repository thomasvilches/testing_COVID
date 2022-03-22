#module covid19abm
using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED  LAT2 PRE2 ASYMP2 MILD2 MISO2 INF2 IISO2 HOS2 ICU2 REC2 DED2 LAT3 PRE3 ASYMP3 MILD3 MISO3 INF3 IISO3 HOS3 ICU3 REC3 DED3 LAT4 PRE4 ASYMP4 MILD4 MISO4 INF4 IISO4 HOS4 ICU4 REC4 DED4 LAT5 PRE5 ASYMP5 MILD5 MISO5 INF5 IISO5 HOS5 ICU5 REC5 DED5 LAT6 PRE6 ASYMP6 MILD6 MISO6 INF6 IISO6 HOS6 ICU6 REC6 DED6 UNDEF
Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    health_status::HEALTH = SUS
    swap::HEALTH = UNDEF
    swap_status::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    wentTo::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    vac_status::Int8 = 0 ##

    got_inf::Bool = false
    herd_im::Bool = false
    hospicu::Int8 = -1
    ag_new::Int16 = -1
    hcw::Bool = false
    days_vac::Int64 = -1
    first_one::Bool = false
    strain::Int16 = -1
    index_day::Int64 = 1
    relaxed::Bool = false
    recovered::Bool = false
    vaccine::Symbol = :none
    vaccine_n::Int16 = 0
    protected::Int64 = 0
    days_recovered::Int64 = -1
    boosted::Bool = false
    n_boosted::Int64 = 0
    recvac::Int64 = 0 # 1 - rec , 2 - vac ... this field shows which immunity will be used for protection

    vac_eff_inf::Array{Array{Array{Float64,1},1},1} = [[[0.0]]]
    vac_eff_symp::Array{Array{Array{Float64,1},1},1} = [[[0.0]]]
    vac_eff_sev::Array{Array{Array{Float64,1},1},1} = [[[0.0]]]

    workplace_idx::Int64 = -1
    school_idx::Int64 = -1

    #### for testing

    daysisolation::Int64 = 999
    days_after_detection::Int64 = 999
    positive::Bool = false
    days_for_pcr::Int64 = -1
    daysinf::Int64 = 999
    tookpcr::Bool = false
    nra::Int64 = 0
    pcrprob::Float64 = 0.0
    test::Bool = false
    
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :ontario
    calibration::Bool = false
    calibration2::Bool = false 
    start_several_inf::Bool = true
    modeltime::Int64 = 435
    initialinf::Int64 = 20
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 1.0  ## percent of people practice self-isolation
    τinf::Int64 = 0
    fsevere::Float64 = 1.0 #
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vac_com_dec_max::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    vac_com_dec_min::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    herd::Int8 = 0 #typemax(Int32) ~ millions
    file_index::Int16 = 0
    nstrains::Int16 = 6
    strain::Int64 = 1 ## which strain is spreading
    strain_test::Int64 = 1 #position of the tuple
    
    #the cap for coverage should be 90% for 65+; 95% for HCW; 80% for 50-64; 60% for 16-49; and then 50% for 12-15 (starting from June 1).
    #comor_comp::Float64 = 0.7 #prop comorbidade tomam
    vaccinating::Bool = true #vaccinating?
    
    ##Alpha - B.1.1.7
    sec_strain_trans::Float64 = 1.5#1.5 #transmissibility of second strain
    
    ## Gamma - P.1
    third_strain_trans::Float64 = 1.6 #transmissibility of third strain
    
    ## Delta - B.1.617.2
    fourth_strain_trans::Float64 = 1.3 #transmissibility compared to second strain strain

    ## Iota - B.1.526
    fifth_strain_trans::Float64 = 1.35 #transmissibility of fifth strain

    ##OMICRON
    rel_trans_sixth::Float64 = 1.35
    sixth_strain_trans::Float64 = rel_trans_sixth*sec_strain_trans*fourth_strain_trans #transmissibility of sixth strain
    reduction_sev_omicron::Float64 = 0.752 ##reduction of severity compared to Delta

    mortality_inc::Float64 = 1.3 #The mortality increase when infected by strain 2

    vaccine_proportion::Vector{Float64} = [0.78; 0.22; 0.0]
    vaccine_proportion_2::Vector{Float64} = [0.78; 0.22; 0.0]
    vac_period::Array{Int64,1} = [21;28;999]
    
    #=------------ Vaccine Efficacy ----------------------------=#
    booster_after::Array{Int64,1} = [180;180;999]
    n_boosts::Int64 = 1
    min_age_booster::Int64 = 16
    #=------------ Vaccine Efficacy ----------------------------=#
    days_to_protection::Array{Array{Array{Int64,1},1},1} = [[[14;21],[0;7]],[[14;21],[0;7]],[[14]]]
    vac_efficacy_inf::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.46;0.46],[0.46;0.861]],[[0.295;0.295],[0.295;0.895]],[[0.368;0.368],[0.368;0.75]],[[0.416;0.416],[0.416;0.85]],[[0.46;0.46],[0.46;0.861]],[[0.16;0.16],[0.16;0.33]]],#booster efficacy  for omicron changed in vac_time function
    [[[0.843;0.843],[0.843,0.964]],[[0.901;0.901],[0.901,0.984]],[[0.67;0.67],[0.67;0.77]],[[0.77;0.77],[0.77,0.867]],[[0.843;0.843],[0.843,0.964]],[[0.16;0.16],[0.16,0.428]]],
    [[[0.61]],[[0.56]],[[0.488]],[[0.496]],[[0.61]],[[0.488]]]]#### 50:5:80

    vac_efficacy_symp::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.63;0.65],[0.65;0.93]],[[0.67;0.7],[0.7;0.89]],[[0.64;0.68],[0.68;0.82]],[[0.57;0.59],[0.59;0.92]],[[0.63;0.65],[0.65;0.93]],[[0.616;0.616],[0.616;0.69]]], #booster efficacy  for omicronmchanged in vac_time function
    [[[0.63;0.7],[0.7,0.96]],[[0.82;0.83],[0.83,0.92]],[[0.75;0.74],[0.74;0.89]],[[0.7;0.69],[0.69,0.95]],[[0.63;0.7],[0.7,0.96]],[[0.678;0.678],[0.678,0.69]]], #### 50:5:80
    [[[0.921]],[[0.88]],[[0.332]],[[0.68]],[[0.921]],[[0.332]]]] #### 50:5:80
    
    vac_efficacy_sev::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.77;0.88],[0.88;0.98]],[[0.82;0.87],[0.87;0.96]],[[0.84;0.87],[0.87;0.96]],[[0.81;0.81],[0.81;0.97]],[[0.77;0.88],[0.88;0.98]],[[0.676;0.676],[0.676;0.81]]], #booster efficacy for omicron changed in vac_time function
    [[[0.66;0.7],[0.7,0.97]],[[0.8;0.82],[0.82,0.95]],[[0.88;0.95],[0.95;0.95]],[[0.9;0.91],[0.91,0.98]],[[0.66;0.7],[0.7,0.97]],[[0.744;0.744],[0.744,0.81]]],#### 50:5:80
    [[[0.921]],[[0.816]],[[0.34]],[[0.781]],[[0.921]],[[0.34]]]]#### 50:5:80

    # ----- Recovery efficacy ----- #
    #https://www.nejm.org/doi/full/10.1056/NEJMc2200133
    # using infection the same as symptoms
    rec_eff_inf::Vector{Float64} = [1.0;0.902;0.857;0.92;1.0;0.56]
    rec_eff_symp::Vector{Float64} = [1.0;0.902;0.857;0.92;1.0;0.56]
    rec_eff_sev::Vector{Float64} = [1.0;0.694;0.88;1.0;1.0;0.878]
    
    time_change_contact::Array{Int64,1} = [1;map(y-> 95+y,0:3);map(y->134+y,0:9);map(y->166+y,0:13);map(y->199+y,0:35)]
    change_rate_values::Array{Float64,1} = [1.0;map(y-> 1.0-0.01*y,1:4);map(y-> 0.96-(0.055/10)*y,1:10);map(y-> 0.90+(0.1/14)*y,1:14);map(y-> 1.0-(0.34/36)*y,1:36)]
    contact_change_rate::Float64 = 1.0 #the rate that receives the value of change_rate_values
    contact_change_2::Float64 = 1.0 ##baseline number that multiplies the contact rate

    relaxed::Bool = false
    relaxing_time::Int64 = 215 ### relax measures for vaccinated
    status_relax::Int16 = 2
    relax_after::Int64 = 1

    turnon::Int64 = 0

    day_inital_vac::Int64 = 104 ###this must match to the matrices in matrice code
    time_vac_kids::Int64 = 253
    time_vac_kids2::Int64 = 428
    using_jj::Bool = false

    #one waning rate for each efficacy? For each strain? I can change this structure based on that
    reduce_days::Int64 = 0
    waning::Int64 = 1
    ### after calibration, how much do we want to increase the contact rate... in this case, to reach 70%
    ### 0.5*0.95 = 0.475, so we want to multiply this by 1.473684211
    proportion_contacts_workplace::Float64 = 0.0
    
    ##for testing
    initial_day_week::Int64 = 1 # 1- Monday ... 7- Sunday
    testing_days::Vector{Int64} = [2;4;6]
    days_ex_test::Int64 = 90 ## 3 months without testing
    isolation_days::Int64 = 5 #how many days of isolation after testing
    test_ra::Int64 = 0 #1 - PCR, 2 - Abbott_PanBio 3 - 	BTNX_Rapid_Response	4 - Artron
    scenariotest::Int64 = 0
    size_threshold::Int64 = 100
    extra_booster::Int64 = 0
end

Base.@kwdef mutable struct ct_data_collect
    total_symp_id::Int64 = 0  # total symptomatic identified
    totaltrace::Int64 = 0     # total contacts traced
    totalisolated::Int64 = 0  # total number of people isolated
    iso_sus::Int64 = 0        # total susceptible isolated 
    iso_lat::Int64 = 0        # total latent isolated
    iso_asymp::Int64 = 0      # total asymp isolated
    iso_symp::Int64 = 0       # total symp (mild, inf) isolated
end

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

include("matrices_code.jl")
## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
#const agebraks_vac = @SVector [0:0,1:4,5:14,15:24,25:44,45:64,65:74,75:100]
const BETAS = Array{Float64, 1}(undef, 0) ## to hold betas (whether fixed or seasonal), array will get resized
const ct_data = ct_data_collect()
const waning_factors = waning_factor()
const waning_factors_rec = waning_factor()
export ModelParameters, HEALTH, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix, nra, npcr = main(ip,simnum)            

    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    #ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    range_work = 18:65
    ags = map(x-> x.age in range_work ? 1 : 2,humans)

    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    work = _collectdf(spl[1])
    
    age_groups = [0:4, 5:17, 18:29, 30:39, 40:49, 50:64, 65:74, 75:84, 85:999]
    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    ag7 = _collectdf(spl[7])
    ag8 = _collectdf(spl[8])
    ag9 = _collectdf(spl[9])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum); insertcols!(ag8, 1, :sim => simnum); 
    insertcols!(ag9, 1, :sim => simnum); insertcols!(work, 1, :sim => simnum);
    

    pos = findall(y-> y in (11,22,33,44,55,66),hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,100)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7, work = work,
    vector_dead=vector_ded,nra=nra,npcr=npcr)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    
    v_strain = [1;4;6]
    p.strain_test = findfirst(x-> x == p.strain,v_strain)

    vac_rate_1::Matrix{Int64} = vaccination_rate_1(sim)
    vac_rate_2::Matrix{Int64} = vaccination_rate_2(sim)
    vac_rate_booster::Vector{Int64} = booster_doses()
    
    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 


    herd_immu_dist_4(sim,1)
    distribute_vaccine(vac_rate_1,vac_rate_2,vac_rate_booster)
    insert_infected(PRE, p.initialinf, 4, ip.strain)[1]
    
    ## save the preisolation isolation parameters
   

    # split population in agegroups 
    grps = get_ag_dist()
    workplaces = create_workplace()
    #schools = create_schools()
    
    initial_dw::Int64 = p.initial_day_week

    nra::Vector{Int64} = zeros(Int64,p.modeltime)
    npcr::Vector{Int64} = zeros(Int64,p.modeltime)

    testing_group::Vector{Int64} = select_testing_group(workplaces)
  

    
    # start the time loop
    for st = 1:p.modeltime
        #= 
        if length(p.time_change_contact) >= count_change && p.time_change_contact[count_change] == st ###change contact pattern throughout the time
            setfield!(p, :contact_change_rate, p.change_rate_values[count_change])
            count_change += 1
        end
        =#
        nra[st],npcr[st] = testing(testing_group,initial_dw)
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,workplaces,schools,initial_dw,sim)
       
        sw = time_update() ###update the system

        initial_dw += 1
        if initial_dw > 7
            initial_dw = 1
        end

        # end of day
    end
    
    
    return hmatrix, nra, npcr## return the model state as well as the age groups. 
end
export main

##### creating workplaces
function work_size() #https://www150.statcan.gc.ca/t1/tbl1/en/cv.action?pid=3310039501
    breaks = [1:4,5:9,10:19,20:49,50:99,100:199,200:499,500:1000]
    #s = [748387, 233347, 152655, 99732, 32889, 14492, 7119, 2803]
    #s/=sum(s)
   # s = [0.5795052593106524,0.18068968828208243,0.11820672374061501,0.07722637956240554,0.025467236167207672,0.01122172113883589,0.00551251951334341,0.0021704722848576454]
    aux = Distributions.Categorical(@SVector [0.5795052593106524,0.18068968828208243,0.11820672374061501,0.07722637956240554,0.025467236167207672,0.01122172113883589,0.00551251951334341,0.0021704722848576454])

    return aux,breaks
end

### I can change it to split basic, mid, and high
function create_workplace() 
    pos = findall(x-> x.age in 20:65,humans)
    N = length(pos)
    
    probs,breaks = work_size()
    vaux = map(y-> rand(breaks[rand(probs)]),1:Int(round(p.popsize/10.0)))
    vvaux = cumsum(vaux)
    aux = findfirst(y-> y > N, vvaux)

    aux == nothing && error("increase the number of workplaces")

    vaux = vaux[1:aux]
    vvaux = vvaux[1:aux]

    vaux[end] = N-vvaux[end-1]
    

    samples::Vector{Vector{Int64}} = map(y-> [y],1:length(vaux))

    for i = 1:length(samples)
        xx = sample(pos,vaux[i],replace=false)
        pos = setdiff(pos,xx)
    
        for j in xx
            humans[j].workplace_idx = i
        end
        samples[i] = deepcopy(xx)
    end


    return samples
end 

function select_testing_group(workplaces::Vector{Vector{Int64}})
    ### we can change this function to implement different test strategies
    if p.scenariotest == 1
        
        wpr = findall(x-> length(x) >= p.size_threshold,workplaces)
        grp = vcat(workplaces[wpr]...)
    elseif p.scenariotest == 2
        wpr = findall(x-> length(x) >= p.size_threshold,workplaces)
        grp = vcat(workplaces[wpr]...)
    elseif p.scenariotest == 3
        grp = 1:length(humans)
    elseif p.scenariotest == 4
        grp = 1:length(humans)
    elseif p.scenariotest == 0
        grp = []
    else
        error("error in testing group")
    end
    
    for i in grp
        humans[i].test = true
    end

    return grp
end

function testing(grp,dayweek)
    npcr::Int64 = 0
    nra::Int64 = 0

    if p.scenariotest == 1
        for x in humans[grp]
            x.days_for_pcr -= 1
            if x.isovia == :symp
                
                if !x.tookpcr#x.daysisolation == 0
                    x.days_for_pcr = rand(1:2)
                    npcr+=1
                    x.pcrprob = _get_prob_test(x,1)
                    x.tookpcr = true
                elseif x.days_for_pcr == 0
                    if rand() > x.pcrprob
                        x.daysisolation = 999
                        x.tookpcr = false
                    end
                end
            end
        end
    elseif p.scenariotest == 2
        for x in humans
            x.days_for_pcr -= 1
            if x.isovia == :symp
                
                if !x.tookpcr#x.daysisolation == 0
                    x.days_for_pcr = rand(1:2)
                    npcr+=1
                    x.pcrprob = _get_prob_test(x,1)
                    x.tookpcr = true
                elseif x.days_for_pcr == 0
                    if rand() > x.pcrprob
                        x.daysisolation = 999
                        x.tookpcr = false
                        x.days_after_detection = 999
                    end
                end
            elseif x.test
                if dayweek in p.testing_days && x.days_after_detection > p.days_ex_test
                    pp = _get_prob_test(x,p.test_ra)
                    if rand() < pp
                        
                        x.positive = true
                        _set_isolation(x,true,:test)
                        
                        x.tookpcr = true
                        x.days_for_pcr = rand(1:2)
                        npcr+=1
                        x.pcrprob = _get_prob_test(x,1)
                        
                    end
                    x.nra += 1
                    nra += 1
                elseif x.days_for_pcr == 0
                    if rand() > x.pcrprob
                        x.daysisolation = 999
                        x.tookpcr = false
                        x.days_after_detection = 999
                    end
                end
            end
        end
    elseif p.scenariotest == 3
        for x in humans
            x.days_for_pcr -= 1
            if x.isovia == :symp
                
                if !x.tookpcr#x.daysisolation == 0
                    x.days_for_pcr = rand(1:2)
                    npcr+=1
                    x.pcrprob = _get_prob_test(x,1)
                    x.tookpcr = true
                elseif x.days_for_pcr == 0
                    if rand() > x.pcrprob
                        x.daysisolation = 999
                        x.tookpcr = false
                    end
                end
            
            end
        end
    elseif p.scenariotest == 4
        for x in humans
            x.days_for_pcr -= 1
            if x.isovia == :symp
                
                pp = _get_prob_test(x,p.test_ra)
                if rand() < pp
                    x.positive = true
                    _set_isolation(x,true,:test)
                end
            
            end
        end
    elseif p.scenariotest == 0

    else
        error("no scenario")
    end


    return nra,npcr
end



function waning_immunity(x::Human)
    index = Int(floor(x.days_vac/7))
    if index > 0
        if index <= size(waning_factors,1)
            waning = [waning_factors[index,x.vaccine_n]^p.waning; waning_factors[index,x.vaccine_n+2]^p.waning]
        else
            waning = [waning_factors[end,x.vaccine_n]^p.waning; waning_factors[end,x.vaccine_n+2]^p.waning]
        end
    else
        waning = [1.0;1.0]
    end

    return waning
end

function vac_update(x::Human)
    
    if x.vac_status == 1
        #x.index_day == 2 && error("saiu com indice 2")
        if x.days_vac == p.days_to_protection[x.vaccine_n][x.vac_status][1]#14
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vaccine_n][x.vac_status]),x.index_day+1)
        elseif x.days_vac == p.days_to_protection[x.vaccine_n][x.vac_status][x.index_day]#14
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vaccine_n][x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.waning = waning_immunity(x)
        x.days_vac += 1

    elseif x.vac_status == 2
        if x.days_vac == p.days_to_protection[x.vaccine_n][x.vac_status][1]#0
            x.protected = 1
            x.index_day = min(length(p.days_to_protection[x.vaccine_n][x.vac_status]),x.index_day+1)

        elseif x.days_vac == p.days_to_protection[x.vaccine_n][x.vac_status][x.index_day]#7
            x.protected = x.index_day
            x.index_day = min(length(p.days_to_protection[x.vaccine_n][x.vac_status]),x.index_day+1)
        end
        if !x.relaxed
            x.relaxed = p.relaxed &&  x.vac_status >= p.status_relax && x.days_vac >= p.relax_after ? true : false
        end
        x.waning = waning_immunity(x)
        x.days_vac += 1
    end
   
end
function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    # resize and update the BETAS constant array
    #init_betas()

    # resize the human array to change population size
    resize!(humans, p.popsize)
end
export reset_params, reset_params_default


## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findall(x-> r[x] == inth && r[x] != r[x-1],2:length(r))
        idx = idx .+ 1
        #idx = findfirst(x -> x == inth, r)
        if idx !== nothing
            for i in idx 
                timevec[i] += 1
            end
        end
    end
    return timevec
end


function herd_immu_dist_4(sim::Int64,strain::Int64)
    rng = MersenneTwister(200*sim)
    vec_n = zeros(Int32,6)
    N::Int64 = 0
    if p.herd == 5
        vec_n = [9; 148; 262;  68; 4; 9]
        N = 5

    elseif p.herd == 10
        vec_n = [32; 279; 489; 143; 24; 33]

        N = 9

    elseif p.herd == 20
        vec_n = [71; 531; 962; 302; 57; 77]

        N = 14
    elseif p.herd == 30
        vec_n = [105; 757; 1448; 481; 87; 122]

        N = 16
    elseif p.herd == 50
        vec_n = map(y->y*5,[32; 279; 489; 143; 24; 33])

        N = 16
    elseif p.herd == 0
        vec_n = [0;0;0;0;0;0]
       
    else
        vec_n = map(y->Int(round(y*p.herd/10)),[32; 279; 489; 143; 24; 33])
        N = 16
    end

    vprob::Vector{Float64} = vector_probs()
    dprob = Distributions.Categorical(vprob)
    for g = 1:6
        pos = findall(y->y.ag_new == g && y.health == SUS,humans)
        n_dist = min(length(pos),Int(floor(vec_n[g]*p.popsize/10000)))
        pos2 = sample(rng,pos,n_dist,replace=false)
        for i = pos2
            humans[i].strain = strain
            humans[i].swap = strain == 1 ? REC : REC2
            humans[i].swap_status = REC
            move_to_recovered(humans[i])
            r = rand()
            day = rand(dprob)
            humans[i].days_recovered = day
            humans[i].sickfrom = INF
            humans[i].herd_im = true
        end
    end
    return N
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin
       #= :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])
        :newyorkcity   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])=#
        :ontario => Distributions.Categorical(@SVector [0.04807822, 0.10498712, 0.12470340, 0.14498051, 0.13137129, 0.12679091, 0.13804896, 0.10292032, 0.05484776, 0.02327152])
        :canada => Distributions.Categorical(@SVector [0.04922255,0.10812899,0.11792442,0.13956709,0.13534216,0.12589012,0.13876094,0.10687438,0.05550450,0.02278485])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag

function comorbidity(ag::Int16)

    a = [4;19;49;64;79;999]
    g = findfirst(x->x>=ag,a)
    prob = [0.05; 0.1; 0.28; 0.55; 0.74; 0.81]

    com = rand() < prob[g] ? 1 : 0

    return com    
end
export comorbidity


function initialize() 
    agedist = get_province_ag(p.prov)
    agebraksnew = [0:4,5:14,15:24,25:34,35:44,45:54,55:64,65:74,75:84,85:99]
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        agn = rand(agedist)
        x.age = rand(agebraksnew[agn]) 
        x.ag = findfirst(y-> x.age in y, agebraks)
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g
        x.exp = 999  ## susceptible people don't expire.
        #x.dur = sample_epi_durations() # sample epi periods   
        x.comorbidity = comorbidity(x.age)
        # initialize the next day counts (this is important in initialization since dyntrans runs first)
        get_nextday_counts(x)
        
    end
end
export initialize

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag,strain) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
    aux_pre = [PRE;PRE2;PRE3;PRE4;PRE5;PRE6]
    aux_lat = [LAT;LAT2;LAT3;LAT4;LAT5;LAT6]
    aux_mild = [MILD;MILD2;MILD3;MILD4;MILD5;MILD6]
    aux_inf = [INF;INF2;INF3;INF4;INF5;INF6]
    aux_asymp = [ASYMP;ASYMP2;ASYMP3;ASYMP4;ASYMP5;ASYMP6]
    aux_rec = [REC;REC2;REC3;REC4;REC5;REC6]
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.strain = strain
            x.first_one = true
            x.dur = sample_epi_durations(x)
            if x.strain > 0
                if health == PRE
                    x.swap = aux_pre[x.strain]
                    x.swap_status = PRE
                    move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
                elseif health == LAT
                    x.swap = aux_lat[x.strain]
                    x.swap_status = LAT
                    move_to_latent(x)
                elseif health == MILD
                    x.swap =  aux_mild[x.strain] 
                    x.swap_status = MILD
                    move_to_mild(x)
                elseif health == INF
                    x.swap = aux_inf[x.strain]
                    x.swap_status = INF
                    move_to_infsimple(x)
                elseif health == ASYMP
                    x.swap = aux_asymp[x.strain]
                    x.swap_status = ASYMP
                    move_to_asymp(x)
                elseif health == REC 
                    x.swap_status = REC
                    x.swap = aux_rec[x.strain]
                    move_to_recovered(x)
                else 
                    error("can not insert human of health $(health)")
                end
            else
                error("no strain insert inf")
            end
            
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function time_update()
    # counters to calculate incidence

    lat_v = zeros(Int64,p.nstrains)
    pre_v = zeros(Int64,p.nstrains)
    asymp_v = zeros(Int64,p.nstrains)
    mild_v = zeros(Int64,p.nstrains)
    miso_v = zeros(Int64,p.nstrains)
    inf_v = zeros(Int64,p.nstrains)
    infiso_v = zeros(Int64,p.nstrains)
    hos_v = zeros(Int64,p.nstrains)
    icu_v = zeros(Int64,p.nstrains)
    rec_v = zeros(Int64,p.nstrains)
    ded_v = zeros(Int64,p.nstrains)

    unvac_r::Int64 = 0
    unvac_nr::Int64 = 0
    vac_1::Int64 = 0
    vac_2::Int64 = 0
    vac_3::Int64 = 0
    
    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
         
        x.daysinf += 1
        x.days_recovered += 1 # We don't care about this up to the recovery
        x.days_after_detection += 1 #we don't care about this untill the individual is detected
        x.daysisolation += 1
        if x.iso && x.daysisolation >= p.isolation_days && !(x.health_status in (HOS,ICU,DED))
            _set_isolation(x,false,:null)
            x.positive = false
        end
        if x.tis >= x.exp             
            @match Symbol(x.swap_status) begin
                :LAT  => begin 
                    move_to_latent(x); 
                    lat_v[x.strain] += 1; 
                    if x.vac_status == 1
                        vac_1 += 1
                    elseif x.vac_status == 2
                        if x.boosted
                            vac_3 += 1
                        else
                            vac_2 += 1
                        end
                    else
                        if x.recovered
                            unvac_r += 1
                        else
                            unvac_nr += 1
                        end
                    end
                end
                :PRE  => begin move_to_pre(x); pre_v[x.strain] += 1; end
                :ASYMP => begin move_to_asymp(x); asymp_v[x.strain] += 1; end
                :MILD => begin move_to_mild(x); mild_v[x.strain] += 1; end
                :MISO => begin move_to_miso(x); miso_v[x.strain] += 1; end
                :INF  => begin move_to_inf(x); inf_v[x.strain] +=1; end    
                :IISO => begin move_to_iiso(x); infiso_v[x.strain] += 1; end
                :HOS  => begin move_to_hospicu(x); hos_v[x.strain] += 1; end 
                :ICU  => begin move_to_hospicu(x); icu_v[x.strain] += 1; end
                :REC  => begin move_to_recovered(x); rec_v[x.strain] += 1; end
                :DED  => begin move_to_dead(x); ded_v[x.strain] += 1; end
                _    => begin dump(x); error("swap expired, but no swap set."); end
            end
        end
        # run covid-19 functions for other integrated dynamics. 
        #ct_dynamics(x)
        # get the meet counts for the next day 
        get_nextday_counts(x)
        if p.vaccinating
            vac_update(x)
        end
       
    end

    #= 
        (lat,lat2,lat3,lat4,lat5,lat6) = lat_v
        (mild,mild2,mild3,mild4,mild5,mild6) = mild_v
        (miso,miso2,miso3,miso4,miso5,miso6) = miso_v
        (inf,inf2,inf3,inf4,inf5,inf6) = inf_v
        (infiso,infiso2,infiso3,infiso4,infiso5,infiso6) = infiso_v
        (hos,hos2,hos3,hos4,hos5,hos6) = hos_v
        (icu,icu2,icu3,icu4,icu5,icu6) = icu_v
        (rec,rec2,rec3,rec4,rec5,rec6) = rec_v
        (ded,ded2,ded3,ded4,ded5,ded6) = ded_v
    =#
    #return (lat, mild, miso, inf, infiso, hos, icu, rec, ded,lat2, mild2, miso2, inf2, infiso2, hos2, icu2, rec2, ded2,lat3, mild3, miso3, inf3, infiso3, hos3, icu3, rec3, ded3, lat4, mild4, miso4, inf4, infiso4, hos4, icu4, rec4, ded4, lat5, mild5, miso5, inf5, infiso5, hos5, icu5, rec5, ded5, lat6, mild6, miso6, inf6, infiso6, hos6, icu6, rec6, ded6)
    return (unvac_r,unvac_nr,vac_1,vac_2,vac_3)
end
export time_update

#@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    if x.isovia == :null
        x.iso = iso 
        x.isovia = via
        x.daysisolation = 0
        x.days_after_detection = 0
    elseif !iso
        x.iso = iso 
        x.isovia = via
    end

end
function sample_epi_durations(y::Human)
    # when a person is sick, samples the 
    
    if y.strain == 4
        lat_dist = Distributions.truncated(LogNormal(1.249, 0.649), 3.5, 7) # truncated between 3.5 and 7
        pre_dist = Distributions.truncated(Gamma(1.015, 1.975), 0.8, 2.2)#truncated between 0.8 and 2.2
    elseif y.strain == 6
        lat_dist = Distributions.truncated(LogNormal(0.99, 0.64), 3, 7) # truncated between 3 and 7
        pre_dist = Distributions.truncated(Gamma(1.015, 1.975), 0.8, 2.2)#truncated between 0.8 and2.2

    else
        lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
        pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    end

    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    aux = y.vac_status > 1 || y.recovered ? p.reduce_days : 0

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = max(Int.(ceil.(rand(asy_dist)))-aux,1)
    infs = max(Int.(ceil.(rand(inf_dist)))-aux,1)

    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = x.swap
    x.health_status = x.swap_status
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
    age_thres = [4, 19, 49, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)

    
    aux_red = 0.0

    if x.recovered
        index = Int(floor(x.days_recovered/7))
        if x.recvac == 1

            if index > 0
                if index <= size(waning_factors_rec,1)
                    aux = waning_factors_rec[index,3]#*(1-aux_red)
                else
                    aux = waning_factors_rec[end,3]#*(1-aux_red)
                end
            else
                aux = 1.0#*(1-aux_red)
            end
            aux = p.rec_eff_symp[x.strain]*aux
        elseif x.recvac == 2

            if x.vac_status*x.protected > 0

                aux = x.vac_eff_symp[x.strain][end][end]*x.waning[2]
                
            else
                if index > 0
                    if index <= size(waning_factors_rec,1)
                        aux = waning_factors_rec[index,3]
                    else
                        aux = waning_factors_rec[end,3]
                    end
                else
                    aux = 1.0
                end
                aux = p.rec_eff_symp[x.strain]*aux
            end
        else
            error("move to latent recvac")
        end
    else
        aux = x.vac_status*x.protected > 0 ? x.vac_eff_symp[x.strain][x.vac_status][x.protected]*x.waning[2] : 0.0
    end
    auxiliar = (1-aux)
 
    if rand() < (symp_pcts[g])*auxiliar

        aux_v = [PRE;PRE2;PRE3;PRE4;PRE5;PRE6]
        x.swap = aux_v[x.strain]
        x.swap_status = PRE
        
    else
        aux_v = [ASYMP;ASYMP2;ASYMP3;ASYMP4;ASYMP5;ASYMP6]
        x.swap = aux_v[x.strain]
        x.swap_status = ASYMP
        
    end
    x.wentTo = x.swap
    x.got_inf = true
    ## in calibration mode, latent people never become infectious.
    if p.calibration && !x.first_one
        x.swap = LAT 
        x.exp = 999
    end 
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = x.swap  
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
   
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    if x.strain == 1 || x.strain == 3 || x.strain == 5 
        θ = (0.95, 0.9, 0.85, 0.6, 0.2)  # percentage of sick individuals going to mild infection stage
    elseif x.strain == 2 || x.strain == 4 || x.strain == 6
        θ = (0.89, 0.78, 0.67, 0.48, 0.04)
        #= if x.strain == 4
            θ = map(y-> max(0,1-(1-y)*1.88),θ)
        end =#
    else
        error("no strain in move to pre")
    end  # percentage of sick individuals going to mild infection stage
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period
    ##########

    aux_red = 0.0
    if x.recovered
        index = Int(floor(x.days_recovered/7))
        aux_red = 0.0#x.strain == 6 ? p.reduction_sev_omicron : 0.0

        if x.recvac == 1

            if index > 0
                if index <= size(waning_factors_rec,1)
                    aux = waning_factors_rec[index,3]#*(1-aux_red)
                else
                    aux = waning_factors_rec[end,3]#*(1-aux_red)
                end
            else
                aux = 1.0#*(1-aux_red)
            end
            aux = p.rec_eff_sev[x.strain]*aux
        elseif x.recvac == 2

            if x.vac_status*x.protected > 0
                aux = x.vac_eff_sev[x.strain][end][end]*x.waning[2]
            else
                if index > 0
                    if index <= size(waning_factors_rec,1)
                        aux = waning_factors_rec[index,3]
                    else
                        aux = waning_factors_rec[end,3]
                    end
                else
                    aux = 1.0
                end
                aux = p.rec_eff_sev[x.strain]*aux
            end
        end
    else
        if x.vac_status*x.protected > 0
            
            aux_red = 0.0#x.strain == 6 ? p.reduction_sev_omicron : 0.0
            aux = x.vac_eff_sev[x.strain][x.vac_status][x.protected]*x.waning[2]
        else
            aux = 0.0
            aux_red = 0.0
        end

    end
    auxiliar = (1-aux)*(1-aux_red)
    
    if rand() < (1-θ[x.ag])*auxiliar
        aux_v = [INF;INF2;INF3;INF4;INF5;INF6]
        x.swap = aux_v[x.strain]
        x.swap_status = INF
    else 
        aux_v = [MILD;MILD2;MILD3;MILD4;MILD5;MILD6]
        x.swap = aux_v[x.strain]
        x.swap_status = MILD
    end
    # calculate whether person is isolated
    #rand() < p.fpreiso && _set_isolation(x, true, :pi)
end
export move_to_pre

function move_to_mild(x::Human)
    ## transfers human h to the mild infection stage for γ days
   
    x.health = x.swap 
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[4]
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    
    #x.swap = x.strain == 1 ? REC : REC2
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, staying in MILD is same as MISO since contacts will be limited. 
    # we still need the separation of MILD, MISO because if x.iso is false, then here we have to determine 
    # how many days as full contacts before self-isolation
    # NOTE: if need to count non-isolated mild people, this is overestimate as isolated people should really be in MISO all the time
    #   and not go through the mild compartment 
    
    if x.iso
        aux_v = [MISO;MISO2;MISO3;MISO4;MISO5;MISO6]
        x.swap = aux_v[x.strain]
        x.swap_status = MISO
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
        
    elseif rand() < p.fmild
        aux_v = [MISO;MISO2;MISO3;MISO4;MISO5;MISO6]
        x.swap = aux_v[x.strain]
        x.swap_status = MISO
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
    end
end
export move_to_mild

function move_to_miso(x::Human)
    ## transfers human h to the mild isolated infection stage for γ days
    x.health = x.swap
    x.health_status = x.swap_status
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0 
    x.exp = x.dur[4] - p.τmild  ## since tau amount of days was already spent as infectious
    
    _set_isolation(x, true, :symp)
   
end
export move_to_miso

function move_to_infsimple(x::Human)
    ## transfers human h to the severe infection stage for γ days 
    ## simplified function for calibration/general purposes
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[4]
    aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
    x.swap = aux_v[x.strain]
    x.swap_status = REC
    #x.swap = x.strain == 1 ? REC : REC2
    _set_isolation(x, false, :null) 
end

function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    # h = prob of hospital, c = prob of icu AFTER hospital    
    comh = 0.98
    if x.strain == 1 || x.strain == 3 || x.strain == 5
        h = x.comorbidity == 1 ? comh : 0.04 #0.376
        c = x.comorbidity == 1 ? 0.396 : 0.25

    elseif x.strain == 2 || x.strain == 4 || x.strain == 6
        if x.age <  20
            h = x.comorbidity == 1 ? comh : 0.05*1.07*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.07 : 0.25*1.07
        elseif x.age >= 20 && x.age < 30
            h = x.comorbidity == 1 ? comh : 0.05*1.29*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.29 : 0.25*1.29
        elseif  x.age >= 30 && x.age < 40
            h = x.comorbidity == 1 ? comh : 0.05*1.45*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.45 : 0.25*1.45
        elseif  x.age >= 40 && x.age < 50
            h = x.comorbidity == 1 ? comh : 0.05*1.61*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.61 : 0.25*1.61
        elseif  x.age >= 50 && x.age < 60
            h = x.comorbidity == 1 ? comh : 0.05*1.58*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.58 : 0.25*1.58
        elseif  x.age >= 60 && x.age < 70
            h = x.comorbidity == 1 ? comh : 0.05*1.65*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.65 : 0.25*1.65
        elseif  x.age >= 70 && x.age < 80
            h = x.comorbidity == 1 ? comh : 0.05*1.45*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.45 : 0.25*1.45
        else
            h = x.comorbidity == 1 ? comh : 0.05*1.60*1 #0.376
            c = x.comorbidity == 1 ? 0.396*1.60 : 0.25*1.60
        end
        #= 
        H. Scribner, Omicron variant leads to less severe symptoms, deaths, new study says. Deseret News (2021) (January 6, 2022).
        UK Health Security Agency, “Technical briefing: Update on hospitalisation and vaccine effectiveness for Omicron VOC-21NOV-01 (B.1.1.529)” (2021).
        C. M aslo, et al., Characteristics and Outcomes of Hospitalized Patients in South Africa During the COVID-19 Omicron Wave Compared With Previous Waves. JAMA (2021) https:/doi.org/10.1001/jama.2021.24868.
        =#
        if x.strain == 4
            if !x.recovered && x.vac_status < 2
                h = h*2.65 #2.26 #https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00475-8/fulltext
            elseif x.recovered || x.boosted #for booster, it is an assumption
                h = h/p.hosp_red #
            end
        elseif x.strain == 6

            h = h*(1-0.3*p.reduction_sev_omicron) # 0.7
            c = c*(1-0.381)#
            if x.recovered || x.boosted #for booster, it is an assumption
                h = h/p.hosp_red
            end
        elseif x.strain == 2
            if x.recovered || x.boosted #for booster, it is an assumption
                h = h/p.hosp_red 
            end
        else
            error("in hospitalization")
        end
    else
        error("no strain in movetoinf")
    end
    
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)

    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
   

    time_to_hospital = Int(round(rand(Uniform(2, 5)))) # duration symptom onset to hospitalization
   	
    x.health = x.swap
    x.health_status = x.swap_status
    x.swap = UNDEF
    
    x.tis = 0 
    if rand() < h     # going to hospital or ICU but will spend delta time transmissing the disease with full contacts 
        x.exp = time_to_hospital
        if rand() < c
            aux_v = [ICU;ICU2;ICU3;ICU4;ICU5;ICU6]
            x.swap = aux_v[x.strain]
            x.swap_status = ICU
            #x.swap = x.strain == 1 ? ICU : ICU2
        else
            aux_v = [HOS;HOS2;HOS3;HOS4;HOS5;HOS6]
            x.swap = aux_v[x.strain]
            x.swap_status = HOS
            #x.swap = x.strain == 1 ? HOS : HOS2
        end
       
    else ## no hospital for this lucky (but severe) individual 
        aux = 0.0#(p.mortality_inc^Int(x.strain==2 || x.strain == 4))
        #aux = x.strain == 4 || x.strain == 6 ? aux*0.0 : aux
        if x.iso || rand() < p.fsevere 
            x.exp = p.τsevere  ## 1 day isolation for severe cases 
            aux_v = [IISO;IISO2;IISO3;IISO4;IISO5;IISO6]
            x.swap = aux_v[x.strain]
            x.swap_status = IISO
            #x.swap = x.strain == 1 ? IISO : IISO2
        else
            if rand() < mh[gg]*aux
                x.exp = x.dur[4] 
                aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
                x.swap = aux_v[x.strain]
                x.swap_status = DED
            else 
                x.exp = x.dur[4]  
                aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
                x.swap = aux_v[x.strain]
                x.swap_status = REC
            end
        end  
    end
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent I -> ?")
end

function move_to_iiso(x::Human)
    ## transfers human h to the sever isolated infection stage for γ days
    x.health = x.swap
    x.health_status = x.swap_status
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)
    
    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for severe cases.
    aux = 0.0#(p.mortality_inc^Int(x.strain==2 || x.strain == 4))
    #aux = x.strain == 4 || x.strain == 6 ? aux*0.0 : aux

    if rand() < mh[gg]*aux
        x.exp = x.dur[4] 
        aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
        x.swap = aux_v[x.strain]
        x.swap_status = DED
    else 
        x.exp = x.dur[4]  
        aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
        x.swap = aux_v[x.strain]
        x.swap_status = REC
    end
    #x.swap = x.strain == 1 ? REC : REC2
    x.tis = 0     ## reset time in state 
    x.exp = x.dur[4] - τinf  ## since 1 day was spent as infectious
    _set_isolation(x, true, :symp)
    
end 

function move_to_hospicu(x::Human)   
    #death prob taken from https://www.cdc.gov/nchs/nvss/vsrr/covid_weekly/index.htm#Comorbidities
    # on May 31th, 2020
    #= age_thres = [24;34;44;54;64;74;84;999]
    g = findfirst(y-> y >= x.age,age_thres) =#
    #https://www.medrxiv.org/content/10.1101/2021.08.24.21262415v1
    aux = [0:4, 5:19, 20:44, 45:54, 55:64, 65:74, 75:84, 85:99]
   
    if x.strain == 1 || x.strain == 3 || x.strain == 5

        mh = [0.001, 0.001, 0.0015, 0.0065, 0.01, 0.02, 0.0735, 0.38]
        mc = [0.002,0.002,0.0022, 0.008, 0.022, 0.04, 0.08, 0.4]

    elseif x.strain == 2  || x.strain == 4  || x.strain == 6
    
        mh = 0.7*[0.0016, 0.0016, 0.0025, 0.0107, 0.02, 0.038, 0.15, 0.66]
        mc = 0.7*[0.0033, 0.0033, 0.0036, 0.0131, 0.022, 0.04, 0.2, 0.70]
        
        if x.strain == 4
            mh = 0.75*mh
            mc = 0.75*mc
        elseif x.strain == 6
            mh = (1-0.9*p.reduction_sev_omicron)*mh
            mc = (1-0.9*p.reduction_sev_omicron)*mc
        end

    else
            error("No strain - hospicu")
    end
    
    gg = findfirst(y-> x.age in y,aux)

    psiH = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17))))
    psiC = Int(round(rand(Distributions.truncated(Gamma(4.5, 2.75), 8, 17)))) + 2
    muH = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15))))
    muC = Int(round(rand(Distributions.truncated(Gamma(5.3, 2.1), 9, 15)))) + 2

    swaphealth = x.swap_status 
    x.health = x.swap ## swap either to HOS or ICU
    x.health_status = x.swap_status
    x.swap = UNDEF
    x.tis = 0
    _set_isolation(x, true, :hosp) # do not set the isovia property here.  

    if swaphealth == HOS
        x.hospicu = 1 
        if rand() < mh[gg] ## person will die in the hospital 
            x.exp = muH 
            aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
            x.swap = aux_v[x.strain]
            x.swap_status = DED
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiH 
            aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
            x.swap = aux_v[x.strain]
            x.swap_status = REC
            #x.swap = x.strain == 1 ? REC : REC2
        end    
    elseif swaphealth == ICU
        x.hospicu = 2 
                
        if rand() < mc[gg] ## person will die in the ICU 
            x.exp = muC
            aux_v = [DED;DED2;DED3;DED4;DED5;DED6]
            x.swap = aux_v[x.strain]
            x.swap_status = DED
            #x.swap = x.strain == 1 ? DED : DED2
        else 
            x.exp = psiC
            aux_v = [REC;REC2;REC3;REC4;REC5;REC6]
            x.swap = aux_v[x.strain]
            x.swap_status = REC
            #x.swap = x.strain == 1 ? REC : REC2
        end
    else
        error("error in hosp")
    end
    
    ## before returning, check if swap is set 
    x.swap == UNDEF && error("agent H -> ?")    
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = h.swap
    h.health_status = h.swap_status
    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    #h.iso = true # a dead person is isolated
    #_set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end

function move_to_recovered(h::Human)
    h.health = h.swap
    h.health_status = h.swap_status
    
    h.recovered = true
    h.days_recovered = 0

    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
    #h.iso = false ## a recovered person has ability to meet others
    h.recvac = 1
    h.daysinf = 999
    
    #_set_isolation(h, false)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function _get_betavalue(sys_time, xhealth) 
    #bf = p.β ## baseline PRE
    #length(BETAS) == 0 && return 0
    bf = p.β#BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.11

    elseif xhealth == MILD || xhealth == MISO 
        bf = bf * 0.44

    elseif xhealth == INF || xhealth == IISO 
        bf = bf * 0.89

    elseif xhealth == ASYMP2
        bf = bf*p.frelasymp*p.sec_strain_trans #0.11

    elseif xhealth == MILD2 || xhealth == MISO2
        bf = bf * 0.44*p.sec_strain_trans

    elseif xhealth == INF2 || xhealth == IISO2 
        bf = bf * 0.89*p.sec_strain_trans

    elseif xhealth == PRE2
        bf = bf*p.sec_strain_trans

    elseif xhealth == ASYMP3
        bf = bf*p.frelasymp*p.third_strain_trans #0.11

    elseif xhealth == MILD3 || xhealth == MISO3
        bf = bf * 0.44*p.third_strain_trans

    elseif xhealth == INF3 || xhealth == IISO3 
        bf = bf * 0.89*p.third_strain_trans

    elseif xhealth == PRE3
        bf = bf*p.third_strain_trans

    elseif xhealth == ASYMP4
        bf = bf*p.frelasymp*p.sec_strain_trans*p.fourth_strain_trans #0.11

    elseif xhealth == MILD4 || xhealth == MISO4
        bf = bf * 0.44*p.sec_strain_trans*p.fourth_strain_trans

    elseif xhealth == INF4 || xhealth == IISO4
        bf = bf * 0.89*p.sec_strain_trans*p.fourth_strain_trans

    elseif xhealth == PRE4
        bf = bf*p.sec_strain_trans*p.fourth_strain_trans
    ############### 5 strain    
    elseif xhealth == ASYMP5
        bf = bf*p.frelasymp*p.fifth_strain_trans #0.11

    elseif xhealth == MILD5 || xhealth == MISO5
        bf = bf * 0.44*p.fifth_strain_trans

    elseif xhealth == INF5 || xhealth == IISO5
        bf = bf * 0.89*p.fifth_strain_trans

    elseif xhealth == PRE5
        bf = bf*p.fifth_strain_trans

    elseif xhealth == ASYMP6
        bf = bf*p.frelasymp*p.sixth_strain_trans #0.11

    elseif xhealth == MILD6 || xhealth == MISO6
        bf = bf * 0.44*p.sixth_strain_trans

    elseif xhealth == INF6 || xhealth == IISO6
        bf = bf * 0.89*p.sixth_strain_trans

    elseif xhealth == PRE6
        bf = bf*p.sixth_strain_trans
    end
    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts
    
    if !x.iso 
        #cnt = rand() < 0.5 ? 0 : rand(1:3)
        aux = x.relaxed ? 1.0*(p.contact_change_rate^p.turnon) : p.contact_change_rate*p.contact_change_2
        cnt = rand(negative_binomials(ag,aux)) ##using the contact average for shelter-in
    else 
        cnt = rand(negative_binomials_shelter(ag,p.contact_change_2))  # expensive operation, try to optimize
    end
    
    if x.health_status == DED
        cnt = 0 
    end
    x.nextday_meetcnt = cnt
    return cnt
end

function dyntrans(sys_time, grps,workplaces,schools,initial_dw,sim)
    totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)
    pos = shuffle(1:length(humans))
    # go through every infectious person
    for x in humans[pos]        
        if x.health_status in (PRE, ASYMP, MILD, MISO, INF, IISO)
            
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts

            if !x.iso && initial_dw ∉ (6,7) ###if no isolated and working days

                if x.workplace_idx > 0 
                    
                    cnts_w = Int(round(cnts*(p.proportion_contacts_workplace)))
                    cnts = cnts-cnts_w 

                    gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                    gpw = [gpw;cnts_w]
        
                    grp_sample = [grps;[workplaces[x.workplace_idx]]]
                    
                else
                    
                    gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                    
                    grp_sample = grps
                    
                end
            else
                gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                grp_sample = grps 
            end


            for (i, g) in enumerate(gpw) 
                meet = rand(grp_sample[i], g)   # sample the people from each group
                # go through each person
                for j in meet 
                    y = humans[j]
                    ycnt = y.nextday_meetcnt    
                    ycnt == 0 && continue

                    y.nextday_meetcnt = y.nextday_meetcnt - 1 # remove a contact
                    totalmet += 1
                    
                    beta = _get_betavalue(sys_time, xhealth)
                    adj_beta = 0 # adjusted beta value by strain and vaccine efficacy
                    if y.health == SUS && y.swap == UNDEF
                        if y.vac_status*y.protected > 0
                            aux = y.vac_eff_inf[x.strain][y.vac_status][y.protected]*y.waning[1]
                            
                        else
                            aux = 0.0
                        end
                         
                        adj_beta = beta*(1-aux)

                    elseif y.health_status == REC && y.swap == UNDEF
                        index = Int(floor(y.days_recovered/7))
                        aux_red = 0.0
                       

                        if y.recvac == 1

                            if index > 0
                                if index <= size(waning_factors_rec,1)
                                    aux = waning_factors_rec[index,1]
                                else
                                    aux = waning_factors_rec[end,1]
                                end
                            else
                                aux = 1.0
                            end
                            #aux = aux*(1-aux_red)
                            aux = p.rec_eff_inf[x.strain]*aux
                        elseif y.recvac == 2

                            if y.vac_status*y.protected > 0
                                aux_vac = y.vac_eff_inf[x.strain][end][end]*y.waning[1]
                                aux = aux_vac*(1-aux_red)
                            else
                                
                                if index > 0
                                    if index <= size(waning_factors_rec,1)
                                        aux = waning_factors_rec[index,1]
                                    else
                                        aux = waning_factors_rec[end,1]
                                    end
                                else
                                    aux = 1.0
                                end
                                
                                aux = p.rec_eff_inf[x.strain]*aux
                            end
                        end

                        adj_beta = beta*(1-aux)
                    end

                    if rand() < adj_beta
                        totalinf += 1
                        y.exp = y.tis   ## force the move to latent in the next time step.
                        y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom
                        y.sickby = x.idx
                        y.strain = x.strain       
                        aux_v = [LAT;LAT2;LAT3;LAT4;LAT5;LAT6]
                        y.swap = aux_v[y.strain]
                        y.swap_status = LAT
                        y.daysinf = 0
                        y.dur = sample_epi_durations(y)
                        #y.swap = y.strain == 1 ? LAT : LAT2
                    end  
                end
            end            
        end
    end
    return totalmet, totalinf
end
export dyntrans

function contact_matrix()
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
     CM[1] = [0.2287, 0.1839, 0.4219, 0.1116, 0.0539]
    CM[2] = [0.0276, 0.5964, 0.2878, 0.0591, 0.0291]
    CM[3] = [0.0376, 0.1454, 0.6253, 0.1423, 0.0494]
    CM[4] = [0.0242, 0.1094, 0.4867, 0.2723, 0.1074]
    CM[5] = [0.0207, 0.1083, 0.4071, 0.2193, 0.2446] 
   
    return CM
end

# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [10.21, 16.793, 13.7950, 11.2669, 8.0027]
    sd = [7.65, 11.7201, 10.5045, 9.5935, 6.9638]
    means = means*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]
end
#const nbs = negative_binomials()
const cm = contact_matrix()
#export negative_binomials, contact_matrix, nbs, cm

export negative_binomials


function negative_binomials_shelter(ag,mult) 
    ## the means/sd here are calculated using _calc_avgag
    means = [2.86, 4.7, 3.86, 3.15, 2.24]
    sd = [2.14, 3.28, 2.94, 2.66, 1.95]
    means = means*mult
    #sd = sd*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]   
end
## internal functions to do intermediate calculations
function _calc_avgag(lb, hb) 
    ## internal function to calculate the mean/sd of the negative binomials
    ## returns a vector of sampled number of contacts between age group lb to age group hb
    dists = _negative_binomials_15ag()[lb:hb]
    totalcon = Vector{Int64}(undef, 0)
    for d in dists 
        append!(totalcon, rand(d, 10000))
    end    
    return totalcon
end
export _calc_avgag

function _negative_binomials_15ag()
    ## negative binomials 15 agegroups
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)
    #0-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70+
    #= AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
     =#
     AgeMean = repeat([14.14],15)#[10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = repeat([10.86],15)#[7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms    
end

#const vaccination_days = days_vac_f()
#const vac_rate_1 = vaccination_rate_1()
#const vac_rate_2 = vaccination_rate_2()
## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
#end # module end
