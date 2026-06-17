#################### R diary for COVID-19 #################################################################################################################################

####################################### Packages used #######################################################################

library(shiny)
library(shinyjs)
library(shinyWidgets)
library(ggplot2)
library(mvtnorm)
library(shinythemes)

###################################### Functions used #######################################################################

########################################################################################################################
#simulating the starting
########################################################################################################################

Sim_Pop_U<-function(N,rad,cluster_pop_radius,area_spread,cluster_inf_radius,loc_cluster,inf_clus) #N:population size
{
  age_l=NULL #lower bound of age group initialisation
  age_u=NULL #upper bound of age group initialisation
  sex=NULL #sex (male/female) initialisation
  working_status=sample(x=c(1:5),N,replace=T,prob = c(18.55 ,6.45 ,22.10 ,20.40 ,32.50)/100) #working status (levels: 1:5) initialisation
  eco_status=sample(x=c(1:5),N,replace=T,prob = c(30.712,45.318,14.981,6.367,2.622)/100) #economic status (levels: 1:5) initialisation
  r_u=rep(2,N) #residence (rural:1, urban:2) initialisation
  abiding_lockdown=rep(0,N) #1:abiding, 0:not abiding
  disease_status=rep(1,N) #1:Susceptible, 2:Presymptomatic, 3:Early Symptomatic, 4:Late Symptomatic, 5:Recovered/Deceased, 6:Asymptomatic
  time_elapsed=rep(0,N) #time after infection occured
  alive=rep(1,N) #1:Alive, 0:Deceased
  test_yn=rep(0,N) #1:Positive 0:Negative
  r=read.table("Urban_age_dist.txt",header=F) #reading urban age dist of West Bengal (Urban) (Source: Census 2011)
  s=sample(c(1:17),N,replace=T,prob=r[,3]) #sampling from the age dist
  p=loc_cluster #locates the coordinates of cluster in the map (slums or overcrowded areas)
  names(p)=c("x","y")
  alloc.prob=0.8/(1-(length(p$x)*pi*cluster_pop_radius^2)/(4*area_spread^2))
  v=c(rep(((1-alloc.prob)/length(p$x)),length(p$x)),alloc.prob) #assigning probability to the clusters (20% of total population is supposed to be in the clusters))
  t=sample(1:length(v),N,replace=T,prob=v) #assigning slum number to an individual
  sex=c(sex,(2-rbinom(N,1,0.529032258)))
  age_l=r[s,1]
  age_u=r[s,2]
  working_status[which(age_l>=5 & age_l<=20)]=sample(c(1:2),sum(age_l>=5 & age_l<=20),replace=T)
  working_status[which(age_l>20)]=sample(c(1:5),sum(age_l>20),replace=T)
  imn_prob = sapply(1:N,function(x){rbinom(1,1,age_sex_immunity(age_l[x],sex[x]))}) #rbinom(N,1,prob = age_sex_immunity(age_l,sex))
  location.matrix<-matrix(0,N,2)
  if(nrow(p) > 0)
  {
    for (i in 1:nrow(p))
    {
      id<-which(t==i)
      location.matrix[id,]<-matrix(as.numeric(p[i,]),length(id),2,byrow=T) + rmvnorm(length(id),c(0,0),sigma=diag(c(cluster_pop_radius/3,cluster_pop_radius/3)))
      eco_status[id]=rbinom(length(id),1,0.5)+1
    }
  }
  id<-which(t==(nrow(p)+1))
  location.matrix[id,]<-matrix(runif(2*length(id),-area_spread,area_spread),length(id),2)
  a=r[s,1]/10
  b=exp(-(a^2)/15)
  trans_23=sapply(b,function(x){sample(c(2:(max(5,floor(8*x)))),1)})
  trans_34=sapply(b,function(x){sample(c(1:(max(3,floor(5*x)))),1)})
  trans_65=sapply(1-b,function(x){sample(c(12:(max(14,floor(21*x)))),1)})
  risk = rbinom(N,1,prob=risk_alive(age_l))
  live_status=1-risk
  trans_45=ifelse(live_status==1,sapply(1-b,function(x){sample(c(4:(max(6,floor(14*x)))),1)}),sapply(b,function(x){sample(c(2:(max(4,floor(7*x)))),1)}))
  dat=data.frame(age_l,age_u,sex,working_status,r_u,location.matrix[,1],location.matrix[,2],abiding_lockdown,disease_status,alive,time_elapsed,trans_23,trans_34,trans_45,trans_65,live_status,eco_status,test_yn,imn_prob)
  plot(dat[,6],dat[,7],pch = 20) #plotting coordinates
  change=NULL
  q=inf_clus #epicenters of initial infection outbreaks
  names(q)=c("x","y")
  if(nrow(q) > 0)
  {
    for(i in 1:min(length(q$x),5))
    {
      x_ind = which(abs(dat[,6] - q$x[i]) < cluster_inf_radius)
      y_ind = which(abs(dat[,7] - q$y[i]) < cluster_inf_radius)
      ind = intersect(x_ind , y_ind) #people inside a neighbourhood of the clicked location (epicenter)
      d_vec = (dat[,6][ind] - rep(q$x[i],length(ind)))^2 + (dat[,7][ind] - rep(q$y[i],length(ind)))^2
      neighbour = ind[which(d_vec < cluster_inf_radius)] #people within 1 unit radius of clicked location (epicenter)
      change=c(change,sample(neighbour,min(floor(0.0004*N),length(neighbour)),replace=F)) #including all peoples inside 1 unit radius of clicked locations
    }
    change = unique(change) #excluding repeatations
    dat[change,9] = 2 #diseased status of those 5% is changed from 1 to 2
    dat[change,11] = 1 #initialising infection time
  }
  dat=as.matrix(dat)
  neighbour_list = lapply(1:N, function(x){setdiff(contact(dat,x,rad),x)}) #list of neighbours of each individual
  ret_list = list("mat"=dat,"neighbours"=neighbour_list) #taking dat matrix and list of neighbours together as a list
  return(ret_list)
}

########################################################################################################################
#gives the mobility of the agents
########################################################################################################################

mobility <- function(dat,index,int_day_mov,area_spread)
{
  final_coor = array(dim = 2)
  f = array(dim = 2)
  init_coor = c(dat[index,6],dat[index,7])
  age_l = dat[index,1]
  sex = dat[index,3]
  working_status = dat[index,4]
  r_u = dat[index,5]
  abiding_lockdown = dat[index,8]
  disease_status = dat[index,9]
  alive = dat[index,10]
  eco_stat = dat[index,17]
  test = dat[index,18]
  if(r_u == 2) #giving movements based on several conditions
  {
    if(disease_status == 1 || disease_status == 2 || disease_status == 6)
    {
      if(age_l == 0 || age_l > 60)
      {
        if(eco_stat < 3 && working_status >= 4)
        {
          final_coor =  rnorm(2,0,0.5/3)
        }
        else if(abiding_lockdown == 1)
        {
          final_coor =  rnorm(2,0,0.01/3)
        }
        else
        {
          final_coor =  rnorm(2,0,0.1/3)
        }
      }
      else if(age_l >= 5 && age_l <= 20)
      {
        final_coor = abiding_lockdown*rnorm(2,0,0.01/3)+(1-abiding_lockdown)*rgamma(2,8,9)*(2*rbinom(2,1,0.5)-1) 
      }
      else
      {
        if(working_status >=4)
        {
          final_coor =  rgamma(2,24,3.5)*(2*rbinom(2,1,0.5)-1)
        }
        else
        {
          if(sex == 1)
          {
            if(abiding_lockdown == 1)
            {
              final_coor =  rnorm(2,0,0.05/3)
            }
            else if(eco_stat < 3)
            {
              final_coor =  rgamma(2,7,18)*(2*rbinom(2,1,0.5)-1)
            }
            else
            {
              final_coor =  rgamma(2,23,5)*(2*rbinom(2,1,0.5)-1)
            }
          }
          else
          {
            if(abiding_lockdown == 1)
            {
              final_coor =  rnorm(2,0,0.01/3)
            }
            else if(eco_stat < 3)
            {
              final_coor =  rnorm(2,0,0.5/3)
            }
            else
            {
              final_coor =  rgamma(2,10,4)*(2*rbinom(2,1,0.5)-1)
            }
          }
        }
      }
    }
    else if(disease_status == 3  || (disease_status == 4 && test == 0))
    {
      if(age_l == 0 || age_l > 60)
      {
        if(eco_stat < 3 && working_status < 4 && abiding_lockdown == 0)
        {
          final_coor =  rnorm(2,0,0.1/3)
        }
        if(eco_stat < 3 && working_status >= 4)
        {
          final_coor =  rnorm(2,0,0.5/3)
        }
        else
        {
          final_coor =  rnorm(2,0,0.01/3)
        }
      }
      else if(age_l >= 5 && age_l <= 20)
      {
        final_coor = rnorm(2,0,ifelse(abiding_lockdown == 1,0.01/3,0.5/3))
      }
      else
      {
        if(working_status >=4)
        {
          final_coor =  rgamma(2,23,5)*(2*rbinom(2,1,0.5)-1)
        }
        else
        {
          if(sex == 1)
          {
            if(abiding_lockdown == 1)
            {
              final_coor =  rnorm(2,0,0.05/3)
            }
            else if(eco_stat < 3)
            {
              final_coor =  rgamma(2,7,18)*(2*rbinom(2,1,0.5)-1)
            }
            else
            {
              final_coor =  rgamma(2,8,9)*(2*rbinom(2,1,0.5)-1)
            }
          }
          else
          {
            if(abiding_lockdown == 1)
            {
              final_coor =  rnorm(2,0,0.01/3)
            }
            else if(eco_stat < 3)
            {
              final_coor =  rnorm(2,0,0.5/3)
            }
            else
            {
              final_coor =  rgamma(2,7,18)*(2*rbinom(2,1,0.5)-1)
            }
          }
        }
      }
    }
    else
    {
      final_coor = c(0,0)
    }
  }
  else
  {
    
  }
  final_coor = (sign(final_coor)*c(min(abs(final_coor[1]),area_spread),min(abs(final_coor[2]),area_spread)))/sqrt(int_day_mov)#if we run this functionn for int_day_mov many times, each time we increment it y X/sqrt(n) times so that
  #the sum of the variables has same mean and variance as that in case of only one movement in a day
  #condition when the final x_co_or is exceeding our 20/20 square, we are reversing the increment to ensure that x-wise it doesn't move out of the box; otherwise final is just the sum of init + incrmnt
  final_coor = (final_coor/5)*area_spread
  f[1]=ifelse((init_coor[1] + final_coor[1] < -area_spread | init_coor[1] + final_coor[1] > area_spread), init_coor[1] - final_coor[1], init_coor[1] + final_coor[1])
  #same for the y_co_or
  f[2]=ifelse((init_coor[2] + final_coor[2] < -area_spread | init_coor[2] + final_coor[2] > area_spread), init_coor[2] - final_coor[2], init_coor[2] + final_coor[2])
  return(f)
}

########################################################################################################################
#returns the index of the neighbours of a particular agent
########################################################################################################################

contact <- function(dat,index,radius)
{
  present_location = dat[index,6:7]
  ind = which(abs(dat[,6] - present_location[1]) < radius & abs(dat[,7] - present_location[2]) < radius) #sorting which of the agents are in the radius of 0.5 units
  neighbour = ind[which((dat[ind,6] - present_location[1])^2 + (dat[ind,7] - present_location[2])^2 < radius^2)] #returns the neighbour of an agent including himself
  return(neighbour)
}

########################################################################################################################
#returns transmission prob for symptomatic neighbours
########################################################################################################################
#

trans_dist<-function(dat,index)  
{
  p=dbeta(seq(0,1,by=1/30),3,7)/3 #we have taken the wts from a +vely skewd distbn as we have assumed in the initial stage of
  x = p[dat[index,11]]            #infection, transmission prob is higher 
  return(x)
}

########################################################################################################################
#returns the disease status of an agent
########################################################################################################################

dis_stat_inf <- function(dat,home_coor,ngbr,index,radius,proportion_asymp,trans_prob_asymp)
{
  r = dat[index,9] #initialising infection status
  im = dat[index,19]
  if(r == 1 & im == 0)
  {
    neighbour = NULL
    if(sum(dat[index,6:7]==home_coor[index,])==2)
    {
      neighbour = ngbr
    }
    else
    {
      neighbour = setdiff(contact(dat,index,radius),index)
    }
    n_dis_stat = dat[neighbour,9] #disease statuses of the people in neighbourhood
    if(sum(n_dis_stat==2 | n_dis_stat==3 | n_dis_stat==4 | n_dis_stat==6) > 0) #there are people (symptomatic) in neighbourhood who can spread
    {
      id.inf = which(n_dis_stat==2 | n_dis_stat==3 | n_dis_stat==4 | n_dis_stat==6)
      inf = n_dis_stat[id.inf]
      pr = rep(trans_prob_asymp,length(inf))
      inf_n = neighbour[which(n_dis_stat==2 | n_dis_stat==3 | n_dis_stat==4)]
      pr[which(inf==2 | inf==3 | inf==4)] = trans_dist(dat,inf_n) #See function "trans_dist"
      eco_change = which(dat[neighbour[id.inf],17] > 3)
      pr[eco_change] = pr[eco_change]/2
      s = sum(rbinom(length(pr),1,pr)) #for each infected people who can infect, a bernoulli rv is generated
      if(s > 0) #if any of the bernoulli rvs are 1 #the person is susceptible #the person will be infected 
      {
        r=4*rbinom(1,1,proportion_asymp)+2 #asymptomatic (stage 6) with 0.25 prob, symptomatic (stage 2) with 0.75 prob
      }
    }
  }
  return(r)
}

########################################################################################################################
#returns the infection time of an agent
########################################################################################################################

time_update<-function(dat,index)
{
  t=dat[index,11]+ifelse((dat[index,9]==1 | dat[index,9]==5), 0, 1) #infection time updation
  return(t)
}

########################################################################################################################
#returns stage transition based on time simulation
########################################################################################################################

dis_stat_stage <- function(dat,index)
{
  curr_stat = dat[index,9] #initialising current status
  if((dat[index,9] == 2 | dat[index,9] == 3 | dat[index,9] == 4) & (dat[index,11] == sum(dat[index,12:(10 + curr_stat)]))) #disease status 2/3/4
  {                                                    #cut-off for moving to status 2/3/4 has been reached
    curr_stat = curr_stat + 1 #status changed to 3/4/5
  }
  else if(dat[index,9] == 6 & dat[index,11] == dat[index,15]) #disease status 6 #cut-off for moving to status 5 reached
  {
    curr_stat = 5 #status changed to 5
  }
  return(curr_stat)
}

########################################################################################################################
#returns immunity probability based on age and sex
########################################################################################################################
#

r=read.table("Urban_age_dist.txt",header=F)
age_given_sex_imo = read.table("age__sex_im0.txt",header = T) # P(age|inf=1,age),column1 = female , column2 = male (probs given in percentage)
age_dist =r[,3] 
sex_prob = 0.529032258
sex2_given_im0 = 3571/13087 # P(sex = female | inf = 1)
j_f_a = (age_given_sex_imo[,1]*(sex2_given_im0/(100*(1-sex_prob))))/age_dist  # (P(age|sex=f,inf=1).P(sex=f|inf=1))/(P(age).P(sex=f))
j_m_a = (age_given_sex_imo[,2]*((1-sex2_given_im0)/(100*sex_prob)))/age_dist  # (P(age|sex=m,inf=1).P(sex=m|inf=1))/(P(age).P(sex=m))
p_im = 1/max(max(j_f_a),max(j_m_a))
imp=cbind((1-(j_m_a*p_im)),(1-(j_f_a*(p_im))))

age_sex_immunity<-function(age,sex)
{
  b=imp #we have defined the immunity probs over the age distribution by using a beta density
  imm_prob = b[((age/5)+1),sex]
  # along with the age factor we have multiplied a sex factor which is 0.7 for females and 0.3 for males as males are less immune as per present data
  return(imm_prob)
}

########################################################################################################################
#returns risk probability of mortality
########################################################################################################################
#
risk_alive<-function(age_l)   #function assigns a indicator of whether a person will live or die if he gets
{
  comorbid = rbinom(length(age_l),1,ifelse(age_l<40, 0.02, 0.35))  #infected based on present conditions of age and comorbidity
  a=exp(2*seq(4,6,by=1/8))/2500000 #among the infected individuals we have generated a risk prob of mortality and while simulating we are drawing a bernoulli with this prob to see if the person will live or die.
  risk_prob=(a[((age_l/5)+1)])*(2*comorbid+1) #For person with comorbidity prob of fatality increases by 3 times over normal chance based on age
  return(risk_prob)
}

########################################################################################################################
#returns the life status of the agent
########################################################################################################################

alive_or_not<-function(dat,index)
{
  live = dat[index,10] #initialising living status
  if(sum(dat[index,9]==4) == length(index)) #disease status 4
  {
    live = ifelse(dat[index,16]==1,1,ifelse(dat[index,11] < rowSums(matrix(dat[index,12:14],length(index))),1,0))
  }
  #col 16 is the indicator of whether a person (in stage 4) will live or not
  #as he will live so the live status will continue to be 1
  #agent will die
  #if the time progression of his infection hasnt reached the cutoff set up by trans_45 col, he will continue to live
  #when he has reached the trans_45 cutoff, his status will change to 0
  return(live)
}


########################################################################################################################
#Function for checking if a person is properly tested
########################################################################################################################

Test_Efficiency<-function(dat,index,sensitivity)
{
  test=rbinom(length(index),1,sensitivity) #test result
  return(test)
}

####################################################################################################################
####################################################################################################################

################################################## APP CODES #######################################################

########################################### UI function ##################################################

ui <- navbarPage("Agent Based Simulation for COVID-19",
                 
  theme = shinytheme("flatly"),
  
  
  
  
  tabPanel("Overview",
           span(p(strong("Introduction")), style="color:#FFFF80; font-size: 12pt"),
           span(p("This is a generic agent based simulator for forecasting disease progression
                  in a small community. Census data were used to assign age, sex, educational status or 
                  occupation to every individual or agent in this community. Information on spread of 
                  COVID-19 infections, chance of detection, transmission potential, course of the disease or 
                  mortality risk obtained from literature were fed into the model. Agents or individuals with 
                  varying attributes behave differently and simulation over time reflects the propagation pattern. 
                  Our simulator tracks every individual, known as 'agent', each day over the entire time period 
                  since the inception of infection till it subsides significantly. Movement or mobility of an agent 
                  depends on age, working status, disease status, intention to follow lockdown etc.
                  We also consider a varying percentage of asymptomatic individuals that play a major role in
                  disease spread. Once an agent dies, s/he is removed from the population and simulation continues 
                  without that agent. We have used comorbidity that depends on age and other factors. 
                  This comorbidity is plugged into the risk probability that indicates the probability that
                  a person at a given stage will live one unit of time."), style="color:#FFFFB3; font-size: 10pt; width: 80%; text-align: justify"),
           br(),
           fluidRow( 
             column( 
           span(p(strong("Choice of parameters")), style="color:#FFFF80; font-size: 12pt"),
           span(p("If you want to run it, first decide on values of some parameters:"),style="color:#FFFFB3; font-size: 10pt; text-align: justify"),
           tags$ul(
           span(tags$li(p(span(strong("Population size (amongst whom you want to forecast):"),style="color:#FFFFB3")," You can vary between a wide range; but for a realistic estimate, it is better to stick to a size of 10000 - 20000.")),
                tags$li(p(span(strong("Proportion of asymptomatic people:"))," In this model, once you fix a value between 0 and 1, it will be treated as constant for the entire simulation period.")),
                tags$li(p(span(strong("Asymptomatic transmission probability:"))," Probability (between 0 and 1) at which one asymptomatic person will infect close contacts.")),
                tags$li(p(span(strong("Radius of infection:"))," It is the radius of a hypothetical circular area, considered as the neighbourhood around the infected person. If an agent enters this neighbourhood, s/he
                will be in close contact of the infected person and is at the risk of getting infected.")),
                tags$li(p(span(strong("Area spread:"))," This is the square shaped space over which the population is distributed. You can fix the area between 10 x 10 and 20 x 20 sq. units.")),
                tags$li(p(strong("Sensitivity:"),"It indicates the conditional probability (between 0 and 1) of identifying an infected person in your setting.")),
                tags$li(p(span(strong("Lock down initiation time:")),"It indicates the commencement of lockdown(s). The input is the starting day of each lockdown. If you prefer to have multiple lockdowns,
                use comma ',' to indicate the day of starting lockdown in phases, from the first day of outbreak.")),
                tags$li(p(span(strong("Common length of lockdown:"))," It is the duration of lockdown (in days). If you want to have multiple lockdowns, in this model, keep duration of lockdown same for different phases.")),
                tags$li(p(span(strong("Population cluster:"))," It is the small pockets in the community with higher population density, usually slums. The number of these clusters can be zero or more, depending on your choice. 
                Fix a reasonable number of clusters between 0 and 5 depending on the area considered for simulation. Population outside these clusters is distributed evenly throughout rest of the area.
                For initiating population cluster, click on the grid and select the option 'Population Cluster'.")),
                tags$li(p(span(strong("Infection cluster:"))," This is the area where first infected case entered the community. Number of infection clusters should be more than one (to get meaningful results) and at most 5. 
                You can choose the infection cluster(s) inside or outside a population cluster.
                For initiating infection cluster, click on the grid and select the option 'Infection Cluster'.")),style="color:#FFFFB3; font-size: 10pt; text-align: justify")),
           width = 6),
           column(
           uiOutput("tb"),width = 6)
           ),
           br(),
           br(),
           span(p("Now you are ready to start the simulation. Click on the 'Run Simulation' tab to get started!!"),style="color:#FFFFB3; font-size: 10pt")
           ),
  
  
  
  tabPanel("Run Simulation",
  
           tags$style(type="text/css",
             ".recalculating {opacity: 1.0;}"),
  
  setBackgroundColor(
    color = "#003333"
    # color = c("#66CCCC", "#003333"),
    # gradient = "radial",
    # direction = c("top", "left")
  ),
  
  useShinyjs(),
  
  sidebarLayout(
    
  
    sidebarPanel(
      radioButtons("start","Let's get started !! (To start click 'Yes')",choices = c("Yes","No"),selected = "No",inline = T),
      # actionButton("reset", "RESET"),
      span(p("Choose the Parameter Values Below:"), style="font-size: 11pt;color:#003333"),
      numericInput("n","Population Size",max = 20000, value = 10000),
      sliderInput("n_sim","Number of Simulations",min = 1,max = 200,value = 150,post = " days"),
      fluidRow(
                column(6,sliderInput("prob_asym","Proportion of Asymptomatic People",min =0.01,max =0.99,value =0.1,step = 0.05)),
                column(6,sliderInput("trans_prob_asym","Asymptomatic Transmission Probability",min =0.01,max =0.99,value =0.2,step = 0.05))
              ),
      fluidRow(
                column(6,sliderInput("rad","Radius of Infection",min = 0.005,max = 0.25,value = 0.15,step = 0.005,post = " units")),
                column(6,sliderInput("area_spread","Area Spread",min =5,max =10,value =10,step = 1))
              ),
      sliderInput("sensitivity","Sensitivity",min = 0.05,max = 1,value = 0.95,step = 0.05),
      radioButtons("lockstart","Do You Want Lockdown Options?",choices = c("Yes","No"),selected = "No",inline = T),
      conditionalPanel("condition = input.lockstart == 'Yes'",
                       span(p("Choose the Following Parameter Values Carefully Based on Simulation Length:"), style="font-size: 11pt;color:#003333"),
                       textInput('vec1', 'Enter Lockdown Initiation Time (comma delimited)', "20,60"),
                       uiOutput("ui"),
                       sliderInput("lockdown_len","Common Length of Lockdowns",min =5,max = 25,value =15,step = 1),
                       uiOutput("ui1")
                       ),
      # radioButtons("start","Let's get started !!",choices = c("Yes","No"),selected = "No",inline = T),
      span(p("Click on Grid for Initiating Clusters. Then Select Appropriate Tab (Population Cluster, Infection Cluster), When Done Press 'START' and Wait (See Progress Bar)",style="font-size: 11pt;color:#003333")),
      actionButton("run", "START"),
      actionButton("reset", "RESET"),
      
      span(p("Make Sure 'Yes' is Selected Above Before Clicking 'START'"), style="font-size: 11pt;color:#003333"),
      width = 4
    ),
    
    
    mainPanel(
      span(textOutput("p1"), style="color:white",align ='center'),
      fluidRow(column(6, align="center", offset = 3, actionButton("action", "Population Cluster"),
                      actionButton("action1", "Infection Cluster"))),
      fluidRow(
               column(12,plotOutput("myplot",click = "plot_click"),style='padding-top:5px; padding-right:10px; padding-left:0px;')
              ),
      
      fluidRow(
               column(4,plotOutput(outputId = "distPlot"),style='padding:0px ;'),
               column(4,plotOutput("spplot"),style='padding:0px ;'),
               column(4,plotOutput("infplot"),style='padding-right:10px;padding-left:0px;')
              )
      )
    ) 
  )
)

############################################# Server unction ##################################################

server <- function(input, output, session) {
  ############################### simulation results construction ##################################
  
  rv = reactiveValues(i = 0)
  c1 = reactiveValues(a1 = "")
  c2 = reactiveValues(a2 = "...")
  c3 = reactiveValues(a3 = "...")
  coors = reactiveValues()
  coors1 = reactiveValues()
  coors$d = data.frame(x = numeric(),y = numeric())
  coors1$d1 = data.frame(x = numeric(),y = numeric())
  
  currentdata <- eventReactive(input$run,{
    if(input$start == "Yes")
    {
      withProgress(message = "App Progress",value = 0,{
        n = input$n #popln size
        n_sim = input$n_sim #no of days of simulation
        cluster_rad = 2 #radius for popln. cluster
        area_sprd = input$area_spread
        cluster_inf_rad = 0.5
        rad = input$rad #radius of infection
        prop_asym = input$prob_asym
        trans_prob_asym = input$trans_prob_asym
        sensitivity=input$sensitivity #sensitivity of the test
        int_day_mov = 1     # input$int_day_mov #no of movements for each person
        lockdown_init_tym = ifelse(input$lockstart =="Yes",sort(as.numeric(unlist(strsplit(input$vec1,",")))),0) 
        no_of_lockdowns = length(lockdown_init_tym)
        lockdown_length = ifelse(input$lockstart =="Yes",rep(input$lockdown_len,no_of_lockdowns),0)
        init = Sim_Pop_U(n,rad,cluster_rad,area_sprd,cluster_inf_rad,coors$d,coors1$d1) #initialisation matrix for n agents
        S = init$mat
        neighourhood = init$neighbours
        start_time <- Sys.time()
        or = S[,6:7]
        
        inf = array(dim=n_sim)
        rec = array(dim=n_sim)
        dec = array(dim=n_sim)
        sus = array(dim=n_sim)
        s_6 = array(dim=n_sim)
        s_2 = array(dim=n_sim)
        s_3 = array(dim=n_sim)
        s_4 = array(dim=n_sim)
        r_0 = array(dim=n_sim)
        tm = array(dim=n_sim)
        plot_list = vector("list", input$n_sim)
        for(g in 1:n_sim)
        {
          plot_list[[g]] = cbind(S[,6],S[,7],S[,9],S[,10])
          if(no_of_lockdowns > 0)
          {
            if(g%in%lockdown_init_tym == TRUE)
            {
              S[,8] = 1
            }
            if(g %in% (lockdown_init_tym + lockdown_length) == TRUE)
            {
              S[,8] = 0
            }
          }
          t_1 = Sys.time()
          change = which(S[,9] != 1 & S[,9] != 5) #agents who are infected at (j-1)-th simulation
          d1_n = which(S[,9] == 4 & S[,18] == 0 & S[,16] == 1) #stage 4 agents in (j-1)-th simulation who are tested negative and have no risk of dying
          c1 = which(S[,9] == 2 | S[,9] == 6) #stage 2 & 6 patients in (j-1)-th simulation
          d = S
          inf_v1 = unique(unlist(lapply(change,function(x){contact(d,x,rad)})))
          susceptible = inf_v1[which(d[inf_v1,9]==1 & d[inf_v1,19]==0)]
          if(length(susceptible) > 0)
          {
            S[susceptible,9] = sapply(susceptible,function(x){dis_stat_inf(d,or,neighourhood[[x]],x,rad,prop_asym,trans_prob_asym)})
          }
          d = S
          t1 = setdiff(which(S[,9] == 2 | S[,9] == 6),c1)
          if(length(t1) > 0)
          {
            S[t1,11] = time_update(d,t1)
          }
          t3 = NULL
          for(i in 1:int_day_mov) #this loop is for interday movement,infection updation and time updation
          {
            d = S
            k1 = which(S[,9] == 2 | S[,9] == 6)
            m = which(d[,9]!= 5 & (d[,9]!=4 | d[,18]!=1))
            if(length(m) > 0)
            {
              S[m,6:7] = t(sapply(m,function(x){mobility(d,x,int_day_mov,area_sprd)}))
            }
            d = S
            inf_v2 = unique(unlist(lapply(which(d[,9]!=1 & d[,9]!=5), function(x){contact(d,x,rad)})))
            susceptible1 = inf_v2[which(d[inf_v2,9]==1 & d[inf_v2,19]==0)]
            if(length(susceptible1)>0)
            {
              S[susceptible1,9] = sapply(susceptible1,function(x){dis_stat_inf(d,or,neighourhood[[x]],x,rad,prop_asym,trans_prob_asym)})
            }
            d = S
            t2 = setdiff(which(S[,9] == 2 | S[,9] == 6),k1)
            t3 = c(t3,t2)
            if(length(t2)>0)
            {
              S[t2,11] = time_update(d,t2)
            }
          }
          S[,6:7] = or   #restoring to original coors
          d = S
          al_n = which(d[,9] == 4)
          if(length(al_n) > 0)
          {
            S[al_n,10] = alive_or_not(d,al_n) #changing the life status of the agents
          }
          d = S
          dss = which(d[,9]!=1 & d[,9]!=5)
          if(length(dss) > 0)
          {
            S[dss,9] = sapply(dss,function(x){dis_stat_stage(d,x)}) #changing the dis_stat column by applying the dis_stat_stage function
          }
          d = S
          l = which(d[,9]==4 & d[,11]==d[,12]+d[,13])
          if(length(l) > 0)
          {
            S[l,18] = rbinom(length(l),1,sensitivity) #Test_Efficiency(d,l,sensitivity) #checking how many stage 4 agents are tested +ve
          }
          d = S
          t4 = setdiff(1:n,union(t1,t3))
          if(length(t4) > 0)
          {
            S[t4,11] = time_update(d,t4)   #changing the time of infection whose time has not been affected yet
          }
          d = S
          test_n = setdiff(which(S[,9] == 4 & S[,18] == 0 & S[,16] == 1),d1_n) #those who are newly not tested in this simulation
          if(length(test_n) > 0)
          {
            S[test_n,16] = 1 - rbinom(length(test_n),1,(0.085+risk_alive(d[test_n,1]))) #we are changing the live_status column of those who initially didnt had any chance
            #of dying but as they are now wrongly tested -ve, so they cant undergo treatment and we again draw a bernoulli with increased prob of dying by  factor 0f 3
          }
          change2 = which(S[,9] != 1 & S[,9] != 5)
          indic = setdiff(change2,change)  #agents who get infected at (j)-th simulation
          r = which(S[,9]==5)
          de = which(S[r,10] == 0)
          s = which(S[,9]==1)
          s2 = which(S[,9]==2)
          s3 = which(S[,9]==3)
          s4 = which(S[,9]==4)
          s6 = which(S[,9]==6)
          inf[g] = length(change2)  #countinng number of infected agents
          r_0[g] = length(indic)/length(change)
          sus[g] = length(s) #countinng number of susceptible agents
          dec[g] = length(de)
          rec[g] = length(setdiff(r,r[de])) #countinr no of stage 5 agents
          s_2[g] = length(s2)
          s_3[g] = length(s3)
          s_4[g] = length(s4)
          s_6[g] = length(s6)
          t_2 = Sys.time()
          tm[g] = t_2 - t_1
          print(paste("g:",g))
          incProgress(1/n_sim)
        }
        
      }
      )
      c2$a2 = "hi"
      data_matr = data.frame(cbind(sus,inf,rec,dec,s_2,s_3,s_4,s_6))
      data_list = list("plot_list" = plot_list,"data" = data_matr,"max" = apply(data_matr,2,max) )
    }
    
  })
  
  ########################################## dynamic plots & texts ################################################
  
  ##################### warning msgs + tutorial video #######################
  
  output$tb <-renderUI({
    h6(
      span(p(strong("Tutorial Video:")), style="color:#FFFF80; font-size: 12pt"),span(p("Here is a short tutorial video on how to use the app:",style="color:#FFFFB3; font-size: 10pt")),
      tags$video(src = "simvid_updated.mp4",width = "700px",height = "450px",type = "video/mp4",controls = "controls")
      )
  })
  
  output$ui <-renderUI({
    t = sort(as.numeric(unlist(strsplit(input$vec1,","))))
    if((t[length(t)] > input$n_sim))
      span(p("Error: Lockdowns Initiation Exceeding No. of Simulations"), style="color:red")
  })
  
  output$ui1 <-renderUI({
    t = sort(as.numeric(unlist(strsplit(input$vec1,","))))
    if((t[length(t)] + input$lockdown_len) > input$n_sim)
      span(p("Error: Lockdown Span Exceeding No. of Simulations"), style="color:red")
  })
  
  ############################### help text #################################
  
  output$p1 <- renderText({
    k = ""
    a = c2$a2
    v = c3$a3
    if(v != k)
    {
      c1$a1 = ifelse(a == "hi","Wait is Over !! Have a Look at the Results (To Restart Press 'RESET')",
                     ifelse(input$start == "Yes",
                            "",k))
    }
    else(c1$a1 = k)
  })
  
  ############################### scatter plot #################################
  
  output$myplot <- renderPlot({
    if(input$start =="Yes")
    {
      par(pty="s")
      plot(coors$d[,1],coors$d[,2],pch = 19, cex = 8,
           col = alpha("red", 0.4),xlim = c(-input$area_spread,input$area_spread),
           ylim =c(-input$area_spread,input$area_spread),
           xlab = "longitude" , ylab = "latitude",
           main = "Population Cluster: Setting Dense Population Clusters in Area \nInfection Cluster: Initiating Cluster Outbreak of Infection (Needed for Non-Null Results)",
           col.main ="grey")
      grid(lty = "solid")
      points(coors1$d1[,1],coors1$d1[,2],pch = 19, cex = 4,
             col = alpha("orange", 0.4),xlim = c(-input$area_spread,input$area_spread),
             ylim =c(-input$area_spread,input$area_spread))
    }
    else{plot.new()}
    if(rv$i > 0) 
    {
      x <- 1:input$n_sim 
      seq_data = currentdata()[[1]]
      S = seq_data[[rv$i]]
      s2 = which(S[,3] == 2)
      s3 = which(S[,3] == 3)
      s4 = which(S[,3] == 4)
      s6 = which(S[,3] == 6)
      r  = which(S[,3] == 5)
      de = which(S[,4][r] == 0)
      s  = which(S[,3] == 1)
      
      S = data.frame(matrix(unlist(S), ncol = 4))
      S[,3][r[de]] = 0
      colnames(S) = c("x","y","Infection","live_status")
      S$Infection = as.factor(S$Infection)
      alpha_vec = rep(0.1,input$n)
      alpha_vec[which(S$Infection != 1)] = 1
      
      sp = ggplot(data = S, mapping = aes(x = x, y = y,color=Infection)) +
        geom_point(alpha = alpha_vec) +
        theme_bw() +  
        theme(aspect.ratio=1) +
        theme(legend.text = element_text(size = 11))+
        theme(legend.title = element_text(size = 12))+
        scale_color_manual(breaks = c("1", "2", "3","4","6","5","0"),values=c("black","blue","orange","violet","yellow","#339900","purple"),
                           labels = c("Susceptible","Pre-symptomatic","Early-symptomatic","Late-symptomatic","Asymptomatic","Recovered","Deceased")) 
      print(sp)
    }
    
  })
  
  ############################### general infection plot #################################
  
  output$distPlot <- renderPlot({
    
    if(rv$i > 0)
    {
      x <- c(1:input$n_sim)[1:rv$i]
      
      l = NULL
      
      if(input$lockstart == "Yes")
      {
        l = sort(as.numeric(unlist(strsplit(input$vec1,","))))
        lock_bar = cbind(l,(l + input$lockdown_len))
      }
      
      
      colors_p = c("Recovered" = "#66CC00", "Active cases" = "red", "Deceased" = "purple")
      
      data_m = currentdata()[[2]]
      data_m = data.frame(matrix(unlist(data_m), ncol = 8))
      data_m = data_m[1:rv$i,]
      colnames(data_m) = c("sus","inf","rec","dec","s_2","s_3","s_4","s_6")
      m = currentdata()[[3]]
      dm = max(m[2],m[3])
      
      
      p = ggplot(data_m) +
        
        geom_line(aes(y=rec, x=x,color = "Recovered"),size = 0.5) +
        
        geom_line(aes(y=inf, x=x,color = "Active cases"),size=0.5) +
        
        geom_line(aes(y=dec, x=x,color = "Deceased"),size=0.5)+
        
        xlim(1,input$n_sim) + 
        ylim(0,dm) +
        labs(title = "Overall Infection Progression",y= "No. of Individuals", x = "No. of Days", color = "Labels") +
        scale_color_manual(values = colors_p) +
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(color = NULL)
      
      if(length(l) > 0)
      {
        if(rv$i > ( lock_bar[1,1]))
        {
          t = which(lock_bar[,1] < rv$i)
          plot_mat = matrix(lock_bar[t,],ncol = 2)
          
          p = p +
            annotate("rect",xmin = plot_mat[t,1], xmax = plot_mat[t,2], ymin = -Inf, ymax = Inf, alpha = 0.1, fill="black") +
            
            annotate(geom="text", x=((plot_mat[t,1]+plot_mat[t,2])/2), y=rep(dm,length(t)), 
                     label = paste0("lockdown: ",1:length(t)),size=4,color="black") 
        }
      }
      
      
      print(p)
    }
    else
    {
      plot.new()
    }
    
  })
  
  ############################### symptomatic progress plot #################################
  
  output$spplot <- renderPlot({
    if(rv$i > 0)
    {
      data_m = currentdata()[[2]]
      
      x <- c(1:input$n_sim)[1:rv$i]
      colors_q = c("Pre-symptomatic" = "blue", "Early-symptomatic" = "orange", "Late-symptomatic" = "6")
      
      l = NULL
      
      if(input$lockstart == "Yes")
      {
        l = sort(as.numeric(unlist(strsplit(input$vec1,","))))
        lock_bar = cbind(l,(l + input$lockdown_len))
      }
      
      data_m = currentdata()[[2]]
      data_m = data.frame(matrix(unlist(data_m), ncol = 8))
      
      data_m = data_m[1:rv$i,]
      colnames(data_m) = c("sus","inf","rec","dec","s_2","s_3","s_4","s_6")
      m = currentdata()[[3]]
      dm = max(m[5],m[7]) + 20
      q = ggplot(data_m) +
        
        geom_line(aes(y=s_2, x=x,color = "Pre-symptomatic"),size=0.5)+
        
        geom_line(aes(y=s_3, x=x,color = "Early-symptomatic"),size=0.5)+
        
        geom_line(aes(y=s_4, x=x,color = "Late-symptomatic"),size=0.5)+
        
        xlim(1,input$n_sim) +
        ylim(0,dm) +
        labs(title = "Infection stages for Symptomatic patients",y= "No. of individuals", x = "No. of days",color = "Labels")+
        scale_color_manual(values = colors_q) +
        labs(color = NULL) +
        theme_classic() +
        theme(legend.position = "bottom")
      
      if(length(l) > 0)
      {
        if(rv$i > ( lock_bar[1,1]  ))
        {
          t = which(lock_bar[,1] < rv$i)
          plot_mat = matrix(lock_bar[t,],ncol = 2)
          
          q = q +
            annotate("rect",xmin = lock_bar[t,1], xmax = lock_bar[t,2], ymin = -Inf, ymax = Inf, alpha = 0.1, fill="black") +
            annotate(geom="text", x=((plot_mat[t,1]+plot_mat[t,2])/2), y=rep(dm,length(t)), 
                     label = paste0("lockdown: ",1:length(t)),size=4,color="black")
        }
      }
      print(q)
      
    }
    else
    {
      plot.new()
    }
    
  })
  
  ############################### symp vs asymp plot #################################
  
  output$infplot <- renderPlot({
    if(rv$i > 0)
    {
      data_m = currentdata()[[2]]
      
      x <- c(1:input$n_sim)[1:rv$i]
      colors_r = c("Symptomatic" = "2", "Asymptomatic" = "yellow")
      
      l = NULL
      
      if(input$lockstart == "Yes")
      {
        l = sort(as.numeric(unlist(strsplit(input$vec1,","))))
        lock_bar = cbind(l,(l + input$lockdown_len))
      }
      
      
      data_m = currentdata()[[2]]
      data_m = data.frame(matrix(unlist(data_m), ncol = 8))
      
      data_m = data_m[1:rv$i,]
      colnames(data_m) = c("sus","inf","rec","dec","s_2","s_3","s_4","s_6")
      m = currentdata()[[3]]
      dm = max(m[2],(m[2]-m[8]))
      r = ggplot(data_m) +
        
        geom_line(aes(y=s_6, x=x,color = "Asymptomatic"),size=0.5)+
        
        geom_line(aes(y=(inf - s_6), x=x,color = "Symptomatic"),size=0.5)+
        xlim(1,input$n_sim) +
        ylim(0,dm) +
        labs(title = "Symptomatic vs Asymptomatic",y = "No. of individuals", x = "No. of days",color = "Labels")+
        scale_color_manual(values = colors_r) +
        theme_classic() +
        labs(color = NULL) +
        theme(legend.position = "bottom")
      
      if(length(l) > 0)
      {
        if(rv$i > ( lock_bar[1,1]  ))
        {
          t = which(lock_bar[,1] < rv$i)
          plot_mat = matrix(lock_bar[t,],ncol = 2)
          
          r = r +
            annotate("rect",xmin = lock_bar[t,1], xmax = lock_bar[t,2], ymin = -Inf, ymax = Inf, alpha = 0.1, fill="black") +
            annotate(geom="text", x=((plot_mat[t,1]+plot_mat[t,2])/2), y=rep(dm,length(t)), 
                     label = paste0("lockdown: ",1:length(t)),size=4,color="black")
        }
      }
      
      print(r)
    }
    else
    {
      plot.new()
    } 
  })
  
  #################################### execution part ########################################
  
  observeEvent(input$run,{ 
    rv$i = 0
    if(input$start == "Yes")
    {
      observe(  
        {
          isolate({ rv$i = rv$i + 1 })
          
          if (isolate(rv$i) < input$n_sim )
          {
            invalidateLater(10, session)
          }
        })}
  })
  
  observeEvent(input$reset, {
    session$reload()
    return()
  })
  
  observe({
    shinyjs::hide("action1")
    shinyjs::hide("action")
    if(input$start == "Yes")
    {show("action")
      show("action1")}
  })
  
  observeEvent(input$run,{
    shinyjs::hide("action")
    shinyjs::hide("action1")
  })
  
  observeEvent(input$action, {
    add_row = data.frame(x = input$plot_click$x,y = input$plot_click$y)
    coors$d = rbind(coors$d, add_row)
  })
  
  observeEvent(input$action1, {
    add_row1 = data.frame(x = input$plot_click$x,y = input$plot_click$y)
    coors1$d1 = rbind(coors1$d1, add_row1)
  })
}

############################################## Run the application ######################################################

shinyApp(ui = ui, server = server)


#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
