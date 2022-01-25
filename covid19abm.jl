module covid19abm
using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
include("matrices_code.jl")
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

    daysafterpositive::Int64 = 999
    positive::Bool = false
    days_for_pcr::Int64 = -1
    daysinf::Int64 = -1
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345       
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :canada
    calibration::Bool = false
    calibration2::Bool = false 
    start_several_inf::Bool = true
    modeltime::Int64 = 435
    initialinf::Int64 = 20
    initialhi::Int64 = 20 ## initial herd immunity, inserts number of REC individuals
    τmild::Int64 = 0 ## days before they self-isolate for mild cases
    fmild::Float64 = 1.0  ## percent of people practice self-isolation
    τinf::Int64 = 0
    fsevere::Float64 = 1.0 #
    eldq::Float64 = 0.0 ## complete isolation of elderly
    eldqag::Int8 = 5 ## default age group, if quarantined(isolated) is ag 5.  
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    fctcapture::Float16 = 0.0 ## how many symptomatic people identified
    #vaccine_ef::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vac_com_dec_max::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    vac_com_dec_min::Float16 = 0.0 # how much the comorbidity decreases the vac eff
    herd::Int8 = 0 #typemax(Int32) ~ millions
    file_index::Int16 = 0
    nstrains::Int16 = 6
    
    #the cap for coverage should be 90% for 65+; 95% for HCW; 80% for 50-64; 60% for 16-49; and then 50% for 12-15 (starting from June 1).
    #comor_comp::Float64 = 0.7 #prop comorbidade tomam

    
    vaccinating::Bool = true #vaccinating?
   
    red_risk_perc::Float64 = 1.0 #relative isolation in vaccinated individuals
    days_Rt::Array{Int64,1} = [100;200;300] #days to get Rt

    ##Alpha - B.1.1.7
    sec_strain_trans::Float64 = 1.5#1.5 #transmissibility of second strain
    ins_sec_strain::Bool = true #insert second strain?
    initialinf2::Int64 = 1 #number of initial infected of second strain
    time_sec_strain::Int64 = 125 #when will the second strain introduced -- Jan 3

    ## Gamma - P.1
    ins_third_strain::Bool = true #insert third strain?
    initialinf3::Int64 = 5 #number of initial infected of third strain
    time_third_strain::Int64 = 999 #when will the third strain introduced - P1 March 20
    third_strain_trans::Float64 = 1.6 #transmissibility of third strain
    
    ## Delta - B.1.617.2
    ins_fourth_strain::Bool = true #insert fourth strain?
    initialinf4::Int64 = 1 #number of initial infected of fourth strain
    time_fourth_strain::Int64 = 173 #when will the fourth strain introduced
    fourth_strain_trans::Float64 = 1.3 #transmissibility compared to second strain strain

    ## Iota - B.1.526
    ins_fifth_strain::Bool = true #insert fifth strain?
    initialinf5::Int64 = 1 #number of initial infected of fifth strain
    time_fifth_strain::Int64 = 7 #when will the fifth strain introduced
    fifth_strain_trans::Float64 = 1.35 #transmissibility of fifth strain

    ##OMICRON
    ins_sixth_strain::Bool = true #insert third strain?
    initialinf6::Int64 = 1 #number of initial infected of sixth strain
    time_sixth_strain::Int64 = 999 #when will the sixth strain introduced
    rel_trans_sixth::Float64 = 1.0
    sixth_strain_trans::Float64 = rel_trans_sixth*sec_strain_trans*fourth_strain_trans #transmissibility of sixth strain
    reduction_sev_omicron::Float64 = 0.7 ##reduction of severity compared to Delta

    mortality_inc::Float64 = 1.3 #The mortality increase when infected by strain 2

    vaccine_proportion::Vector{Float64} = [0.59;0.33;0.08]
    vaccine_proportion_2::Vector{Float64} = [0.63;0.37;0.0]
    vac_period::Array{Int64,1} = [21;28;999]
    booster_after::Array{Int64,1} = [180;180;999]
    booster_after_bkup::Array{Int64,1} = [120;120;999]
    change_booster_eligibility::Int64 = 999
    n_boosts::Int64 = 1
    time_first_to_booster::Int64 = 9999
    min_age_booster::Int64 = 16
    reduction_omicron::Float64 = 0.0
    reduction_reduction::Float64 = 0.0
    #=------------ Vaccine Efficacy ----------------------------=#
    
    #=------------ Vaccine Efficacy ----------------------------=#
    days_to_protection::Array{Array{Array{Int64,1},1},1} = [[[14],[0;7]],[[14],[0;14]],[[14]]]
    vac_efficacy_inf::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.46],[0.6;0.861]],[[0.295],[0.6;0.895]],[[0.368],[0.48;0.736]],[[0.368],[0.48;0.64]],[[0.46],[0.6;0.861]],[[0.368],[0.48;0.64]]],
    [[[0.61],[0.61,0.935]],[[0.56],[0.56,0.86]],[[0.488],[0.488;0.745]],[[0.496],[0.496,0.76]],[[0.61],[0.61,0.935]],[[0.496],[0.496,0.76]]],
    [[[0.61]],[[0.56]],[[0.488]],[[0.496]],[[0.61]],[[0.488]]]]#### 50:5:80

    vac_efficacy_symp::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.57],[0.66;0.94]],[[0.536],[0.62;0.937]],[[0.332],[0.66;0.94]],[[0.335],[0.62;0.88]],[[0.57],[0.66;0.94]],[[0.335],[0.62;0.88]]],
    [[[0.921],[0.921,0.941]],[[0.88],[0.88,0.91]],[[0.332],[0.66;0.94]],[[0.68],[0.68,0.70]],[[0.921],[0.921,0.941]],[[0.68],[0.68,0.70]]], #### 50:5:80
    [[[0.921]],[[0.88]],[[0.332]],[[0.68]],[[0.921]],[[0.332]]]] #### 50:5:80
    
    vac_efficacy_sev::Array{Array{Array{Array{Float64,1},1},1},1} = [[[[0.62],[0.80;0.92]],[[0.541],[0.8;0.94]],[[0.34],[0.68;0.974]],[[0.34],[0.68;0.80]],[[0.62],[0.80;0.92]],[[0.34],[0.68;0.80]]],
    [[[0.921],[0.921,1.0]],[[0.816],[0.816,0.957]],[[0.34],[0.68;0.974]],[[0.781],[0.781,0.916]],[[0.921],[0.921,1.0]],[[0.781],[0.781,0.916]]],#### 50:5:80
    [[[0.921]],[[0.816]],[[0.34]],[[0.781]],[[0.921]],[[0.34]]]]#### 50:5:80


    time_change_contact::Array{Int64,1} = [1;map(y-> 95+y,0:3);map(y->134+y,0:9);map(y->166+y,0:13);map(y->199+y,0:35)]
    change_rate_values::Array{Float64,1} = [1.0;map(y-> 1.0-0.01*y,1:4);map(y-> 0.96-(0.055/10)*y,1:10);map(y-> 0.90+(0.1/14)*y,1:14);map(y-> 1.0-(0.34/36)*y,1:36)]
    contact_change_rate::Float64 = 1.0 #the rate that receives the value of change_rate_values
    contact_change_2::Float64 = 0.5 ##baseline number that multiplies the contact rate

    relaxed::Bool = false
    relaxing_time::Int64 = 215 ### relax measures for vaccinated
    status_relax::Int16 = 2
    relax_after::Int64 = 1

    relax_over::Int64 = 92
    relax_rate::Float64 = (1-contact_change_2)/relax_over
    turnon::Int64 = 0
    time_back_to_normal::Int64 = 999

    day_inital_vac::Int64 = 104 ###this must match to the matrices in matrice code
    time_vac_kids::Int64 = 253
    time_vac_kids2::Int64 = 428
    using_jj::Bool = false

    α::Float64 = 1.0
    α2::Float64 = 0.0
    α3::Float64 = 1.0
    daysofvac::Int64 = 365### check it
    doubledose::Int64 = 999
    doubledose_kids::Int64 = 999
    rate_dd_kids::Float64 = 2.0
    booster_increase::Int64 = 999
    increase_booster_rate::Float64 = 1.0

    scenario::Symbol = :statuscuo

    #one waning rate for each efficacy? For each strain? I can change this structure based on that

    waning::Int64 = 0
    ### after calibration, how much do we want to increase the contact rate... in this case, to reach 70%
    ### 0.5*0.95 = 0.475, so we want to multiply this by 1.473684211
    proportion_contacts_workplace::Float64 = 0.0
    proportion_contacts_school::Float64 = 0.0

    ##for testing
    initial_day_week::Int64 = 1 # 1- Monday ... 7- Sunday
    testing_days::Vector{Int64} = [2;4;6]
    days_ex_test::Int64 = 90 ## 3 months without testing
    isolation_test_days::Int64 = 10 #how many days of isolation after testing
    pcrend::Bool = false #perform a pcr at the end of isolation?
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
    hmatrix, remaining_doses, total_given, unvac_rec,unvac_unrec,vac_1,vac_2,vac_3 = main(ip,simnum)            

    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    #ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    range_work = 18:65
    ags = map(x-> x.age in range_work ? 1 : 2,humans)

    all = _collectdf(hmatrix)
    spl = _splitstate(hmatrix, ags)
    work = _collectdf(spl[1])
    
    age_groups = [0:4, 5:11, 12:17, 18:49, 50:64, 65:79, 80:999]
    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    ag7 = _collectdf(spl[7])
    insertcols!(all, 1, :sim => simnum); insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum); insertcols!(work, 1, :sim => simnum);
    

    coverage1 = length(findall(x-> x.age >= 18 && x.vac_status >= 1,humans))/length(findall(x-> x.age >= 18,humans))
    coverage2 = length(findall(x-> x.age >= 18 && x.vac_status == 2,humans))/length(findall(x-> x.age >= 18,humans))

    coverage12 = length(findall(x-> x.vac_status >= 1,humans))/p.popsize
    coverage22 = length(findall(x-> x.vac_status == 2,humans))/p.popsize

    #### let's count the number of vaccines for each vaccine thaat was given
    aux =  findall(x-> x.vaccine_n == 2, humans)
    n_moderna = length(aux)
    aux =  findall(x-> x.vaccine_n == 1, humans)
    n_pfizer = length(aux)
    aux =  findall(x-> x.vaccine_n == 3, humans)
    n_jensen = length(aux)

    aux =  findall(x-> x.vaccine_n == 2 && x.age in range_work, humans)
    n_moderna_w = length(aux)
    aux =  findall(x-> x.vaccine_n == 1 && x.age in range_work, humans)
    n_pfizer_w = length(aux)
    aux =  findall(x-> x.vaccine_n == 3 && x.age in range_work, humans)
    n_jensen_w = length(aux)

    aux =  findall(x-> x.vaccine_n == 2 && x.vac_status == 2, humans)
    n_moderna_2 = length(aux)
    aux =  findall(x-> x.vaccine_n == 1 && x.vac_status == 2, humans)
    n_pfizer_2 = length(aux)
    aux =  findall(x-> x.vaccine_n == 3 && x.vac_status == 2, humans)
    n_jensen_2 = length(aux)

    aux =  findall(x-> x.vaccine_n == 2 && x.age in range_work && x.vac_status == 2, humans)
    n_moderna_w_2 = length(aux)
    aux =  findall(x-> x.vaccine_n == 1 && x.age in range_work && x.vac_status == 2, humans)
    n_pfizer_w_2 = length(aux)
    aux =  findall(x-> x.vaccine_n == 3 && x.age in range_work && x.vac_status == 2, humans)
    n_jensen_w_2 = length(aux)

    pos = findall(y-> y in (11,22,33,44,55,66),hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,100)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    return (a=all, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7, work = work,
    cov1 = coverage1,cov2 = coverage2,cov12 = coverage12,cov22 = coverage22,vector_dead=vector_ded,
    unvac_nr = unvac_unrec, unvac_r=unvac_rec, vac_1 = vac_1,vac_2 = vac_2, vac_3 = vac_3,
    n_pfizer = n_pfizer, n_moderna = n_moderna, n_jensen = n_jensen, n_pfizer_w = n_pfizer_w, n_moderna_w = n_moderna_w, n_jensen_w = n_jensen_w,
    n_pfizer_2 = n_pfizer_2, n_moderna_2 = n_moderna_2, n_jensen_2 = n_jensen_2, n_pfizer_w_2 = n_pfizer_w_2, n_moderna_w_2 = n_moderna_w_2, n_jensen_w_2 = n_jensen_w_2, remaining = remaining_doses, total_given = total_given)
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
    
    vac_rate_1::Matrix{Int64} = vaccination_rate_1(sim)
    vac_rate_2::Matrix{Int64} = vaccination_rate_2(sim)
    vac_rate_booster::Vector{Int64} = booster_doses()
    vaccination_days::Vector{Int64} = days_vac_f(size(vac_rate_1,1))

    if p.doubledose-p.day_inital_vac < length(vac_rate_booster)
        vac_rate_1[(p.doubledose-p.day_inital_vac):end,:] = 2*vac_rate_1[(p.doubledose-p.day_inital_vac):end,:]
        vac_rate_2[(p.doubledose-p.day_inital_vac):end,:] = 2*vac_rate_2[(p.doubledose-p.day_inital_vac):end,:]
        vac_rate_booster[(p.doubledose-p.day_inital_vac):end] = 2*vac_rate_booster[(p.doubledose-p.day_inital_vac):end]
    end
    
    if p.doubledose_kids-p.day_inital_vac < size(vac_rate_1,1)
        vac_rate_1[(p.doubledose_kids-p.day_inital_vac):end,1:2] = Int.(round.(p.rate_dd_kids*vac_rate_1[(p.doubledose_kids-p.day_inital_vac):end,1:2]))
        vac_rate_2[(p.doubledose_kids-p.day_inital_vac):end,1:2] = Int.(round.(p.rate_dd_kids*vac_rate_2[(p.doubledose_kids-p.day_inital_vac):end,1:2]))
        #vac_rate_booster[p.doubledose_kids:end] = vac_rate_booster[p.doubledose_kids:end]
    end

    if p.booster_increase-p.day_inital_vac < length(vac_rate_booster)
        vac_rate_booster[(p.booster_increase-p.day_inital_vac):end] = Int.(round.(p.increase_booster_rate*vac_rate_booster[(p.booster_increase-p.day_inital_vac):end]))
    end

    
    agebraks_vac::SVector{8, UnitRange{Int64}} = get_breaks_vac()#@SVector [0:0,1:4,5:14,15:24,25:44,45:64,65:74,75:100]

    v_prop,fd_prop,sd_prop = temporal_proportion()
    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 
   
    N = herd_immu_dist_4(sim,1)
    if p.initialinf > 0
        insert_infected(PRE, p.initialinf, 4, 1)[1]
    end
        #findall(x->x.health in (MILD,INF,LAT,PRE,ASYMP),humans)
   
    h_init1 = findall(x->x.health  in (LAT,MILD,MISO,INF,PRE,ASYMP),humans)
    h_init1 = [h_init1]
    h_init2 = []
    h_init3 = []
    h_init4 = []
    ## save the preisolation isolation parameters
    _fpreiso = p.fpreiso
    p.fpreiso = 0

    # split population in agegroups 
    grps = get_ag_dist()
    workplaces = create_workplace()
    schools = create_schools()
    count_change::Int64 = 1
    
    time_vac::Int64 = 1
    time_pos::Int64 = 0
    time_prop::Int64 = 1
    remaining_doses::Int64 = 0
    total_given::Int64 = 0
    count_relax::Int64 = 1

    initial_dw::Int64 = p.initial_day_week

    unvac_rec::Vector{Int64} = zeros(Int64,p.modeltime)
    unvac_unrec::Vector{Int64} = zeros(Int64,p.modeltime)
    vac_1::Vector{Int64} = zeros(Int64,p.modeltime)
    vac_2::Vector{Int64} = zeros(Int64,p.modeltime)
    vac_3::Vector{Int64} = zeros(Int64,p.modeltime)

    testing_group::Int64 = select_testing_group()
  

    if p.vaccinating
        vac_ind::Vector{Vector{Int64}} = vac_selection(sim,16,agebraks_vac)
    else
        time_vac = 9999 #this guarantees that no one will be vaccinated
    end
    # start the time loop
    for st = 1:p.modeltime
        if p.ins_sec_strain && st == p.time_sec_strain ##insert second strain
            insert_infected(PRE, p.initialinf2, 4, 2)[1]
            h_init2 = findall(x->x.health  in (LAT2,MILD2,INF2,PRE2,ASYMP2),humans)
            h_init2 = [h_init2]
        end
        if p.ins_third_strain && st == p.time_third_strain #insert third strain
            insert_infected(PRE, p.initialinf3, 4, 3)[1]
            h_init3 = findall(x->x.health  in (LAT3,MILD3,INF3,PRE3,ASYMP3),humans)
            h_init3 = [h_init3]
        end
        if p.ins_fourth_strain && st == p.time_fourth_strain #insert third strain
            insert_infected(PRE, p.initialinf4, 4, 4)[1]
            h_init4 = findall(x->x.health  in (LAT4,MILD4,INF4,PRE4,ASYMP4),humans)
            h_init4 = [h_init4]
        end
        if p.ins_fifth_strain && st == p.time_fifth_strain #insert third strain
            insert_infected(PRE, p.initialinf5, 4, 5)[1]
        end
        if p.ins_sixth_strain && st == p.time_sixth_strain #insert third strain
            insert_infected(PRE, p.initialinf6, 4, 6)[1]
        end
        if length(p.time_change_contact) >= count_change && p.time_change_contact[count_change] == st ###change contact pattern throughout the time
            setfield!(p, :contact_change_rate, p.change_rate_values[count_change])
            count_change += 1
        end

        if p.vaccinating
            if st == p.time_vac_kids
                vac_ind = vac_selection(sim,12,agebraks_vac)
            elseif st == p.time_vac_kids2
                vac_ind = vac_selection(sim,5,agebraks_vac)
            
            end
        end 

        if st == p.change_booster_eligibility
            p.booster_after = deepcopy(p.booster_after_bkup)
        end
        # start of day
        #println("$st")

        if st == p.relaxing_time ### time that people vaccinated people is allowed to go back to normal
            setfield!(p, :relaxed, true)
        end
        if st >= p.time_back_to_normal && count_relax <= p.relax_over
            #setfield!(p, :contact_change_2, p.contact_change_2+p.relax_rate)
            p.contact_change_2 += p.relax_rate
            count_relax += 1
        end

        if time_pos < length(vaccination_days) && time_vac == vaccination_days[time_pos+1]
            time_pos += 1
        end

        if time_prop < length(v_prop) && st == v_prop[time_prop]
            setfield!(p, :vaccine_proportion, fd_prop[time_prop,:])
            setfield!(p, :vaccine_proportion_2, sd_prop[time_prop,:])
            time_prop += 1
        end

        #= if p.vaccinating && st == p.time_vac_kids 
            vac_ind = vac_selection(sim,12,agebraks_vac)
        end =#

        time_vac += 1
        if time_pos > 0 
            if st >= p.time_first_to_booster
                vac_rate_booster[time_pos+1] += sum(vac_rate_1[time_pos+1,:])
                vac_rate_1[time_pos+1,:] .= 0
            end
            aux_ =  vac_time!(sim,vac_ind,time_pos+1,vac_rate_1,vac_rate_2,vac_rate_booster)
            remaining_doses += aux_[1]
            total_given += aux_[2]
        end

        #auxx = findall(x-> x.vac_status > 0,humans)
        
       #=  if length(auxx) > 0
            aa = sum([humans[x].vac_status for x in auxx]) 
            if aa != total_given
                println("$st $(time_pos+1) erro $total_given $aa")
            end
        end =#

        #println([time_vac length(findall(x-> x.vac_status == 2 && x.age >= 18,humans))])
       
        testing(grp,initial_dw)
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,workplaces,schools,initial_dw,sim)
       
        sw = time_update() ###update the system

        initial_dw += 1
        if initial_dw > 7
            initial_dw = 1
        end

        unvac_unrec[st] = sw[2]
        unvac_rec[st] = sw[1]
        vac_1[st] = sw[3]
        vac_2[st] = sw[4]
        vac_3[st] = sw[5]
        # end of day
    end
    
    
    return hmatrix, remaining_doses, total_given, unvac_rec,unvac_unrec,vac_1,vac_2,vac_3 ## return the model state as well as the age groups. 
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

function select_testing_group()
    ### we can change this function to implement different test strategies
    grp = findall(y-> y.age in 12:65,humans)
    return grp
end

function testing(grp,dayweek)
    npcr::Int64 = 0
    ntest::Int64 = 0
    for i in grp
        x = humans[i]
        if !x.positive 
            #x.daysafterpositive += 1 #need to add one here
            if dayweek in p.testing_days && x.daysafterpositive > p.days_ex_test #check if they will be tested
                test_individual(x)
                ntest+=1
            end
        end
    end

    for x in humans
        
        if x.iso
            if x.daysafterpositive == x.days_for_pcr ## gap between antigen and pcr results
                test_individual_pcr(x, false)
                npcr+=1
            end
            if x.daysafterpositive > p.isolation_test_days
                if p.pcrend #performing a new pcr?
                    test_individual_pcr(x, true)
                    npcr+=1
                else
                    _set_isolation(x, false, :null)
                    x.positive = false
                end
            end
        end
        x.daysafterpositive += 1 #need to add one here
        if x.daysinf >= 0 
            x.daysinf += 1
        end
    end

    return ntest,npcr
end

function test_individual(x::Human)
    pp = _get_prob_test(x)
    if rand() < pp
        x.positive = true
        x.daysafterpositive = 0
        _set_isolation(x, true, :test)
    end
end

function test_individual_pcr(x::Human,ending::Bool)
    pp = _get_prob_test_pcr(x)

    if rand() > pp #the individual is already positive, so, if rand() > pp, we cancel the positivity
        x.positive = false
        x.daysafterpositive = 999
        _set_isolation(x, false, :null)
    else
        x.daysafterpositive = ending ? 0 : x.daysafterpositive
    end
end

##### creating schools
function schools_size()
    size = 100
    return size
end

### I can change it to split basic, mid, and high
function create_schools()
    pos = findall(x-> x.age in 5:18,humans)
    N = length(pos)
    
    vaux = map(y-> schools_size(),1:Int(round(p.popsize/10.0)))
    vvaux = cumsum(vaux)
    aux = findfirst(y-> y > N, vvaux)
    aux == nothing && error("increase the number of schools")

    vaux = vaux[1:aux]
    vvaux = vvaux[1:aux]

    vaux[end] = N-vvaux[end-1]
    

    samples::Vector{Vector{Int64}} = map(y-> [y],1:length(vaux))

    for i = 1:length(samples)
        xx = sample(pos,vaux[i],replace=false)
        pos = setdiff(pos,xx)
    
        for j in xx
            humans[j].school_idx = i
        end
        samples[i] = deepcopy(xx)
    end


    return samples
end 

function waning_immunity(x::Human)
    index = Int(floor(x.days_vac/7))
    if index > 0 && index <= size(waning_factors,1)
        j = x.protected
        for i in 1:length(x.vac_eff_inf)
            x.vac_eff_inf[i][x.vac_status][j] = p.vac_efficacy_inf[x.vaccine_n][i][x.vac_status][j]*(waning_factors[index,x.vaccine_n])^p.waning
            x.vac_eff_symp[i][x.vac_status][j] = p.vac_efficacy_symp[x.vaccine_n][i][x.vac_status][j]*(waning_factors[index,x.vaccine_n+2])^p.waning
            x.vac_eff_sev[i][x.vac_status][j] = p.vac_efficacy_sev[x.vaccine_n][i][x.vac_status][j]*(waning_factors[index,x.vaccine_n+2])^p.waning
        end
    end
end



function vac_selection(sim::Int64,age::Int64,agebraks_vac)
    
    

   
    aux_1 = map(k-> findall(y-> y.age in k && y.age >= age && y.comorbidity == 1,humans),agebraks_vac)
    aux_2 = map(k-> findall(y-> y.age in k && y.age >= age && y.comorbidity == 0,humans),agebraks_vac)

    v = map(x-> [aux_1[x];aux_2[x]],1:length(aux_1))
    
    return v
end


function vac_time!(sim::Int64,vac_ind::Vector{Vector{Int64}},time_pos::Int64,vac_rate_1::Matrix{Int64},vac_rate_2::Matrix{Int64},vac_rate_booster::Vector{Int64})
    aux_states = (MILD, MISO, INF, IISO, HOS, ICU, DED)
    ##first dose
   # rng = MersenneTwister(123*sim)
    ### lets create distribute the number of doses per age group
    
    remaining_doses::Int64 = 0
    total_given::Int64 = 0

    ### Let's calculate the number of available doses per vaccine type
    doses_first::Vector{Int64} = Int.(round.(sum(vac_rate_1[time_pos,:])*(p.vaccine_proportion/sum(p.vaccine_proportion))))
    doses_second::Vector{Int64} = Int.(round.(sum(vac_rate_2[time_pos,:])*(p.vaccine_proportion_2/sum(p.vaccine_proportion_2))))

    if sum(doses_first) < sum(vac_rate_1[time_pos,:])
        r = rand(1:3)
        doses_first[r] += 1
    elseif sum(doses_first) > sum(vac_rate_1[time_pos,:])
        rr=findall(y-> y>0,doses_first)
        r = rand(rr)
        doses_first[r] -= 1
    end

    if sum(doses_second) < sum(vac_rate_2[time_pos,:])
        r = rand(1:2)
        doses_second[r] += 1
    elseif sum(doses_second) > sum(vac_rate_2[time_pos,:])
        rr=findall(y-> y>0,doses_second)
        r = rand(rr)
        doses_first[r] -= 1
    end


    for i in 1:length(vac_ind)
        pos = findall(y-> humans[y].vac_status == 1 && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[i])
        
        l1 = min(vac_rate_2[time_pos,i],length(pos))
        remaining_doses += (vac_rate_2[time_pos,i])

        for j = 1:l1
            if doses_second[1] > 0 && doses_second[2] > 0 
                pos2 = filter(x-> humans[vac_ind[i][x]].vac_status == 1,pos)
            elseif doses_second[1] > 0 
                pos2 = filter(x-> humans[vac_ind[i][x]].vaccine_n == 1 && humans[vac_ind[i][x]].vac_status == 1,pos)
            elseif doses_second[2] > 0 
                pos2 = filter(x-> humans[vac_ind[i][x]].vaccine_n == 2  && humans[vac_ind[i][x]].vac_status == 1,pos)
            else
                error("missing doses")
            end

            if length(pos2) > 0 
                r = rand(pos2)
                x = humans[vac_ind[i][r]]
                x.days_vac = 0
                x.vac_status = 2
                x.index_day = 1
                total_given += 1
                doses_second[x.vaccine_n] -= 1
                remaining_doses -= 1

                if x.recovered
                    index = Int(floor(x.days_recovered/7))

                    if index > 0
                        if index <= size(waning_factors_rec,1)
                            aux = waning_factors_rec[index,1]
                        else
                            aux = waning_factors_rec[end,1]
                        end
                    else
                        aux = 1.0
                    end

                    if aux > p.vac_efficacy_inf[x.vaccine_n][1][x.vac_status][end]
                        x.recvac = 1
                    else
                        x.recvac = 2
                    end
                end


            else
                break
            end
            
        end

        pos = findall(y-> humans[y].vac_status == 0 && !(humans[y].health_status in aux_states),vac_ind[i])
        
        l2 = min(vac_rate_1[time_pos,i],length(pos))
        remaining_doses += (vac_rate_1[time_pos,i])
        for j = 1:l2
            if doses_first[1] > 0
                pos2 = filter(x-> humans[vac_ind[i][x]].vac_status == 0,pos)
            else
                pos2 = filter(x-> humans[vac_ind[i][x]].age >= 18 && humans[vac_ind[i][x]].vac_status == 0,pos)
            end

            if length(pos2) > 0 
                r = rand(pos2)
                x = humans[vac_ind[i][r]]
                x.days_vac = 0
                x.vac_status = 1
                x.index_day = 1
                x.vaccine_n = x.age < 18 ? 1 : sample([1,2,3], Weights(doses_first/sum(doses_first)))
                x.vaccine = [:pfizer;:moderna;:jensen][x.vaccine_n]
                doses_first[x.vaccine_n] -= 1
                remaining_doses -= 1
                total_given += 1

                x.vac_eff_inf = deepcopy(p.vac_efficacy_inf[x.vaccine_n])
                x.vac_eff_symp = deepcopy(p.vac_efficacy_symp[x.vaccine_n])
                x.vac_eff_sev = deepcopy(p.vac_efficacy_sev[x.vaccine_n])


                if x.recovered
                    index = Int(floor(x.days_recovered/7))

                    if index > 0
                        if index <= size(waning_factors_rec,1)
                            aux = waning_factors_rec[index,1]
                        else
                            aux = waning_factors_rec[end,1]
                        end
                    else
                        aux = 1.0
                    end

                    if aux > p.vac_efficacy_inf[x.vaccine_n][1][x.vac_status][end]
                        x.recvac = 1
                    else
                        x.recvac = 2
                    end
                end
            else
                break
            end
            
        end
    end
    ###remaining_doses are given to any individual within the groups that are vaccinated on that day
    v1 = ones(Int64,doses_second[1]+doses_first[1])
    v2 = 2*ones(Int64,doses_second[2]+doses_first[2])
    v3 = 3*ones(Int64,doses_second[3]+doses_first[3])

    v = shuffle([v1;v2;v3])

    for vac in v
        if vac == 1
            pos = map(k->findall(y-> humans[y].vac_status == 1 && humans[y].vaccine_n == vac && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
            pos2 = map(k->findall(y-> humans[y].vac_status == 0 && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
        
        elseif vac == 2
            pos = map(k->findall(y-> humans[y].vac_status == 1 && humans[y].vaccine_n == vac && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
            pos2 = map(k->findall(y-> humans[y].vac_status == 0 && humans[y].age >= 18 && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
        
        else
            pos = map(k->findall(y-> humans[y].vac_status == 1 && humans[y].vaccine_n == vac && humans[y].days_vac >= p.vac_period[humans[y].vaccine_n] && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
            pos2 = map(k->findall(y-> humans[y].vac_status == 0 && humans[y].age >= 18 && !(humans[y].health_status in aux_states),vac_ind[k]),1:length(vac_ind))
        end

        aux = findall(x-> vac_rate_1[time_pos,x] > 0 || vac_rate_2[time_pos,x] > 0, 1:length(vac_ind))
        position = map(k-> vac_ind[k][pos[k]],aux)
        position2 = map(k-> vac_ind[k][pos2[k]],aux)
        r = vcat(position...,position2...)
        
        length(r) == 0 && continue
            
        rr = sample(r)
        x = humans[rr]
        if x.vac_status == 0
            x.days_vac = 0
            x.vac_status = 1
            x.index_day = 1
            x.vaccine_n = vac
            x.vaccine = [:pfizer;:moderna;:jensen][x.vaccine_n]
            remaining_doses -= 1
            total_given += 1
            x.vac_eff_inf = deepcopy(p.vac_efficacy_inf[x.vaccine_n])
            x.vac_eff_symp = deepcopy(p.vac_efficacy_symp[x.vaccine_n])
            x.vac_eff_sev = deepcopy(p.vac_efficacy_sev[x.vaccine_n])
            if x.recovered
                index = Int(floor(x.days_recovered/7))

                if index > 0
                    if index <= size(waning_factors_rec,1)
                        aux = waning_factors_rec[index,1]
                    else
                        aux = waning_factors_rec[end,1]
                    end
                else
                    aux = 1.0
                end

                if aux > p.vac_efficacy_inf[x.vaccine_n][1][x.vac_status][end]
                    x.recvac = 1
                else
                    x.recvac = 2
                end
            end
        elseif x.vac_status == 1
            x.days_vac = 0
            x.vac_status = 2
            x.index_day = 1
            remaining_doses -= 1
            total_given += 1

            if x.recovered
                index = Int(floor(x.days_recovered/7))

                if index > 0
                    if index <= size(waning_factors_rec,1)
                        aux = waning_factors_rec[index,1]
                    else
                        aux = waning_factors_rec[end,1]
                    end
                else
                    aux = 1.0
                end

                if aux > p.vac_efficacy_inf[x.vaccine_n][1][x.vac_status][end]
                    x.recvac = 1
                else
                    x.recvac = 2
                end
            end
            
        else
            error("error in humans vac status - vac time")
        end
    end

   
    ### Let's add booster... those are extra doses, we don't care about missing doses

    pos = findall(y-> y.vac_status == 2 && y.days_vac >= p.booster_after[y.vaccine_n] && y.age >= p.min_age_booster && y.n_boosted < p.n_boosts && !(y.health_status in aux_states),humans)

    l2 = min(vac_rate_booster[time_pos]+remaining_doses,length(pos))
    pos = sample(pos,l2,replace=false)

    for i in pos
        x = humans[i]
        x.days_vac = 0
        x.boosted = true
        x.n_boosted += 1
        x.vac_eff_inf = deepcopy(p.vac_efficacy_inf[x.vaccine_n])
        x.vac_eff_symp = deepcopy(p.vac_efficacy_symp[x.vaccine_n])
        x.vac_eff_sev = deepcopy(p.vac_efficacy_sev[x.vaccine_n])

        if x.recovered
            index = Int(floor(x.days_recovered/7))

            if index > 0
                if index <= size(waning_factors_rec,1)
                    aux = waning_factors_rec[index,1]
                else
                    aux = waning_factors_rec[end,1]
                end
            else
                aux = 1.0
            end

            if aux > p.vac_efficacy_inf[x.vaccine_n][1][x.vac_status][end]
                x.recvac = 1
            else
                x.recvac = 2
            end
        end
    end

   

    return remaining_doses,total_given

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
        waning_immunity(x)
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
        idx = findfirst(x -> x == inth, r)
        if idx !== nothing 
            timevec[idx] += 1
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
            day = findfirst(y-> y > r, vprob)
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
        :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :canada => Distributions.Categorical(@SVector [0.0540, 0.1697, 0.3915, 0.2159, 0.1689])
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
        :ontario => Distributions.Categorical(@SVector [0.0519, 0.1727, 0.3930, 0.2150, 0.1674])
        :newyorkcity   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])
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
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        x.ag = rand(agedist)
        x.age = rand(agebraks[x.ag]) 
        a = [4;19;49;64;79;999]
        g = findfirst(y->y>=x.age,a)
        x.ag_new = g
        x.exp = 999  ## susceptible people don't expire.
        x.dur = sample_epi_durations() # sample epi periods   
        if rand() < p.eldq && x.ag == p.eldqag   ## check if elderly need to be quarantined.
            x.iso = true   
            x.isovia = :qu         
        end
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
        if x.recovered 
            x.days_recovered += 1
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

@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
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
    x.iso = iso 
    x.isovia == :null && (x.isovia = via)

    x.days_for_pcr = x.iso ? rand(1:2) : 999

end

function sample_epi_durations()
    # when a person is sick, samples the 
    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = Int.(ceil.(rand(asy_dist)))
    infs = Int.(ceil.(rand(inf_dist)))
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

    index = Int(floor(x.days_recovered/7))
    aux_red = 0.0

    if x.recovered
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

        elseif x.recvac == 2

            if x.vac_status*x.protected > 0

                aux = x.vac_eff_symp[x.strain][end][end]
                
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
            end
        else
            error("move to latent recvac")
        end
    else
        aux = x.vac_status*x.protected > 0 ? x.vac_eff_symp[x.strain][x.vac_status][x.protected] : 0.0
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

    
    if x.recovered
        index = Int(floor(x.days_recovered/7))
        aux_red = x.strain == 6 ? p.reduction_sev_omicron : 0.0

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

        elseif x.recvac == 2

            if x.vac_status*x.protected > 0
                aux = x.vac_eff_sev[x.strain][end][end]
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
            end
        end
    else
        if x.vac_status*x.protected > 0
            
            aux_red = x.strain == 6 ? p.reduction_sev_omicron : 0.0
            aux = x.vac_eff_sev[x.strain][x.vac_status][x.protected]
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
    aux = x.vac_status > 0 ? p.fmild*p.red_risk_perc : p.fmild

    if x.iso
        aux_v = [MISO;MISO2;MISO3;MISO4;MISO5;MISO6]
        x.swap = aux_v[x.strain]
        x.swap_status = MISO
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
        
    elseif rand() < aux#p.fmild
        aux_v = [MISO;MISO2;MISO3;MISO4;MISO5;MISO6]
        x.swap = aux_v[x.strain]
        x.swap_status = MISO
        #x.swap = x.strain == 1 ? MISO : MISO2  
        x.exp = p.τmild
        x.daysafterpositive = 0
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
    
    _set_isolation(x, true, :mi) 
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
                
        if x.strain == 4
            if !x.recovered && x.vac_status < 2
                h = h*2.26 #https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00475-8/fulltext
            elseif x.recovered || x.boosted #for booster, it is an assumption
                h = h/6.7 #https://www.medrxiv.org/content/10.1101/2021.08.24.21262415v1
            end
        elseif x.strain == 6
            #https://www.deseret.com/coronavirus/2021/12/31/22861222/omicron-variant-less-severe-covid-symptoms-deaths
            #https://jamanetwork.com/journals/jama/fullarticle/2787776?guestAccessKey=919da83d-b6f9-4e05-8de1-05cca4541a59&utm_source=silverchair&utm_medium=email&utm_campaign=article_alert-jama&utm_content=olf&utm_term=123021
            h = h*(1-p.reduction_sev_omicron) # 0.7
            c = c*(1-0.36)#https://jamanetwork.com/journals/jama/fullarticle/2787776?guestAccessKey=919da83d-b6f9-4e05-8de1-05cca4541a59&utm_source=silverchair&utm_medium=email&utm_campaign=article_alert-jama&utm_content=olf&utm_term=123021
            if x.recovered || x.boosted #for booster, it is an assumption
                h = h/6.7 #https://www.medrxiv.org/content/10.1101/2021.08.24.21262415v1
            end
        elseif x.strain == 2
            if x.recovered || x.boosted #for booster, it is an assumption
                h = h/6.7 #https://www.medrxiv.org/content/10.1101/2021.08.24.21262415v1
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
   
    ###prop/(prob de sintoma severo)
    if p.calibration && !p.calibration2
        h =  0#, 0, 0, 0)
        c =  0#, 0, 0, 0)
        mh = (0, 0, 0, 0, 0)
    end

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
        aux = (p.mortality_inc^Int(x.strain==2 || x.strain == 4))
        aux = x.strain == 4 || x.strain == 6 ? aux*0.0 : aux
        if x.iso 
            x.exp = τinf  ## 1 day isolation for severe cases 
            aux_v = [IISO;IISO2;IISO3;IISO4;IISO5;IISO6]
            x.swap = aux_v[x.strain]
            x.swap_status = IISO
        elseif rand() < p.fsevere 
            x.exp = τinf  ## 1 day isolation for severe cases 
            aux_v = [IISO;IISO2;IISO3;IISO4;IISO5;IISO6]
            x.swap = aux_v[x.strain]
            x.swap_status = IISO
            x.daysafterpositive = 0
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
    aux = (p.mortality_inc^Int(x.strain==2 || x.strain == 4))
    aux = x.strain == 4 || x.strain == 6 ? aux*0.0 : aux

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
    x.exp = x.dur[4] - 1  ## since 1 day was spent as infectious
    _set_isolation(x, true, :mi)
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
    
        mh = 0.5*[0.0016, 0.0016, 0.0025, 0.0107, 0.02, 0.038, 0.15, 0.66]
        mc = 0.5*[0.0033, 0.0033, 0.0036, 0.0131, 0.022, 0.04, 0.2, 0.70]
        
        if x.strain == 4
            mh = 1.0*mh
            mc = 1.0*mc
        elseif x.strain == 6
            #https://www.ijidonline.com/action/showPdf?pii=S1201-9712%2821%2901256-X
            # 70% -> reduction
            mh = (1-p.reduction_sev_omicron)*mh
            mc = (1-p.reduction_sev_omicron)*mc
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
    h.daysinf = -1
    
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
                    if x.school_idx > 0
                        ## carefull here
                        cnts_s = Int(round(cnts*(p.proportion_contacts_school)))
                        cnts_w = Int(round(cnts*(p.proportion_contacts_workplace)))
                        cnts = cnts-cnts_w-cnts ##carefull here

                        gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                        gpw = [gpw;cnts_w;cnts_s]
            
                        grp_sample = [grps;[workplaces[x.workplace_idx]];[schools[x.school_idx]]]
                    else
                        cnts_w = Int(round(cnts*(p.proportion_contacts_workplace)))
                        cnts = cnts-cnts_w 

                        gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                        gpw = [gpw;cnts_w]
            
                        grp_sample = [grps;[workplaces[x.workplace_idx]]]
                    end
                else
                    if x.school_idx > 0
                        cnts_s = Int(round(cnts*(p.proportion_contacts_school)))
                        cnts = cnts-cnts_s

                        gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                        gpw = [gpw;cnts_s]
            
                        grp_sample = [grps;[schools[x.school_idx]]]
                    else
                        
                        gpw = Int.(round.(cm[x.ag]*cnts)) # split the counts over age groups
                        
                        grp_sample = grps
                    
                    end
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

                            if x.strain == 6
                                if y.boosted
                                    aux_r = (y.days_vac > 90) ? (1-p.reduction_omicron*(1-p.reduction_reduction)) : 1.0
                                else
                                    aux_r = (y.days_vac > 90 && y.vac_status == 2) ? (1-p.reduction_omicron) : 1.0
                                end
                            else
                                aux_r = 1.0
                            end

                            aux = aux_r*y.vac_eff_inf[x.strain][y.vac_status][y.protected]
                        else
                            aux = 0.0
                        end
                         
                        adj_beta = beta*(1-aux)

                    elseif y.health_status == REC && y.swap == UNDEF
                        index = Int(floor(y.days_recovered/7))

                        if y.vac_status > 0
                            aux_red = (x.strain == 6 && y.days_recovered > 90 && y.days_vac > 90) ? p.reduction_omicron : 0.0
                        else
                            aux_red = (x.strain == 6 && y.days_recovered > 90) ? p.reduction_omicron : 0.0
                        end

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
                            aux_vac = y.vac_status*y.protected > 0 ? (1-p.reduction_reduction) : 1.0
                            aux = aux*(1-aux_red*aux_vac)

                        elseif y.recvac == 2

                            if y.vac_status*y.protected > 0

                                aux_vac = y.vac_eff_inf[x.strain][end][end]
                                aux = aux_vac*(1-aux_red*(1-p.reduction_reduction))
                            else
                                
                                if index > 0
                                    if index <= size(waning_factors_rec,1)
                                        aux = waning_factors_rec[index,1]*(1-aux_red)
                                    else
                                        aux = waning_factors_rec[end,1]*(1-aux_red)
                                    end
                                else
                                    aux = 1.0*(1-aux_red)
                                end
    
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
end # module end
