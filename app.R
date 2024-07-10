library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)

simulatey = function(A1,A2,A3,l,wbar,cbar,alpha,eta,p,q,s,w,r,m,it){
  
  #Independent of t, fixed from the outset
  n=3
  L=1
  B = matrix(rbind(A1,A2,A3),ncol=n)
  #l = c(l1,l2,l3)
  #wbar = c(wbar1,wbar2,wbar3)
  #cbar = c(cbar1,cbar2,cbar3)
  #alpha = c(alpha1,alpha2)
  #eta = c(eta1,eta2,eta3,eta4,eta5,eta6,eta7,eta8)
  
  #Dependent on t
  x = list()
  x$p = p
  x$q = t(q)
  x$s = s
  x$w = w
  x$r = r
  x$m = m
  
  #Finding k's
  kis = x$p*x$s^eta[1:n]
  kw_finder = function(w,L,l,q,eta,n){
    w*(L-sum(l*q))^eta[(2*n+1)]
  }
  k_w = kw_finder(x$w,L,l,x$q,eta,n)
  kr_finder = function(r,m,eta,n){
    r*m[2]^eta[(2*n+2)]
  }
  k_r = kr_finder(x$r,x$m,eta,n)
  k = c(kis,k_w,k_r)
  
  data = c()
  
  #Begin iteration process here after preparing x
  for (j in 1:it){
    #Stats
    x$realwage = alpha[1]*x$m[1]*wbar/sum(x$p*wbar)
    x$capconsu = alpha[2]*x$m[2]*cbar/sum(x$p*cbar)
    
    demand = c()
    for (f in 1:n){
      demand = c(demand,sum(B[f,]*x$q)+x$realwage[f]+x$capconsu[f])
    }
    x$demand = demand
    
    costs = c()
    for (f in 1:n){
      costs = c(costs,(sum(x$p*B[,f])+l[f]*x$w)*x$q[f])
    }
    x$costs = costs
    
    x$profits = x$p*pmin(x$demand,x$q+x$s)-x$costs*(1+x$r)
    x$roi = (x$p*pmin(x$demand,x$q+x$s)-x$costs)/x$costs
    x$employments = l*x$q
    x$agg_employ = sum(x$employments)
    x$agg_demand = sum(x$demand*x$p)
    
    x$h = x$costs/x$q
    C = x$capconsu%*%x$h/sum(x$h*x$q)
    Bsquig = B + C
    x$naturals = l%*%solve(diag(n)-Bsquig)*x$w
    
    #Save state
    data = cbind(data,unlist(x))
    
    #Find dps from (7.7)
    dps = -eta[1:n]*k[1:n]^(-1/eta[1:n])*x$p^(rep(1,n)+1/eta[1:n])*(x$q-x$demand)
    
    #Find dqs from (7.8)
    dqs = eta[n+(1:n)]*x$q*x$profits/(x$costs*(1+x$r))
    
    #Find dmw from (7.1) and (7.11)
    dmw = k[n+1]*x$agg_employ/(L-x$agg_employ)^eta[2*n+1]-alpha[1]*x$m[1]
    
    #Money changes hands
    x$m = x$m + c(dmw,-dmw)
    
    #Capitalist decisions for next period
    x$p = x$p + dps
    x$q = x$q + dqs
    
    #Inventory updates from (7.10)
    x$s = (k[1:n]/x$p)^(1/eta[1:n])
    
    #Wage updates from (7.11)
    x$w = k[n+1]/(L-(sum(l*x$q)))^eta[2*n+1]
    
    #Interest rate updates from (7.12)
    x$r = k[n+2]/x$m[2]^(eta[2*n+2])
    
  }
  #End loop
  
  Y=data.frame(t(data))
  Y
  
}

graffy = function(chunk,names,title="Graph",leg=FALSE){
  data = melt(chunk)
  data$Var2 = factor(data$Var2)
  ggplot(data, aes(x = Var1, y = value, color = Var2, group = Var2)) +
    geom_line(show.legend = leg) +
    ggtitle(title) +
    scale_color_discrete(labels = names) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank())
}

dumbgraffy = function(chunk,title="Graph"){
  ggplot(data.frame(time = 1:length(chunk), value = chunk),
         aes(x = time, y = value)) +
    geom_line()+
    ggtitle(title) +
    theme(axis.title.y=element_blank(),
          axis.title.x=element_blank())
}

graffyrex = function(Y){
  r=dim(Y)[1]
  n=3
  pnames = paste("p",1:n,sep="")
  nnames = paste("naturals",1:n,sep="")
  names = c(pnames,nnames)
  
  pricevals <- Y %>% select(all_of(names))
  
  pricevals = cbind(stack(pricevals[,1:n]),stack(pricevals[,(n+1):(2*n)]))
  names(pricevals) = c("price","fam","TLC","lol")
  
  pricevals %>% 
    ggplot(aes(x=price, y=TLC,group=fam,color=fam)) +
    ggtitle("Prices p vs. Total Labor Costs v\u0303w") +
    geom_point(size=0.6,show.legend=FALSE) +
    geom_function(color = "black",fun=function(a){a}) +
    geom_segment(aes(
      xend=c(price[2:r],NA,price[(r+2):(2*r)],NA,price[(2*r+2):(3*r)],NA), 
      yend=c(TLC[2:r],NA,TLC[(r+2):(2*r)],NA,TLC[(2*r+2):(3*r)],NA)
    ),
    show.legend=FALSE,
    arrow=arrow(length=unit(0.13,"cm"))
    ) 
}


## prepare data structures to create UI programmatically
parms = list(A=matrix(c(0.2,0.2,0,0,0.8,0,0.4,0,0.1),ncol=3,
          dimnames = list(c("Corn needed","Iron needed","Sugar needed"), 
                          c("to make Corn", "to make Iron","to make Sugar"))),
          l=c(0.7,0.6,0.3),
          wbar=c(0.6,0,0.2),
          cbar=c(0.2,0,0.4),
          alpha=c(0.8,0.7),
          eta=c(2,2,2,1,1,1,0.25,2),
          p=c(1,0.8,0.5),q=c(0.01,0.1,0.2),s=c(2,2,3),w=0.5,r=0.03,m=c(0.5,0.5),
          it=40)

server <- function(input, output) {
   
  
  output$everything <- renderPlot({
    
    parms = list()
    parms$A1 = as.numeric(unlist(strsplit(input$A1,",")))
    parms$A2 = as.numeric(unlist(strsplit(input$A2,",")))
    parms$A3 = as.numeric(unlist(strsplit(input$A3,",")))
    parms$l = as.numeric(unlist(strsplit(input$l,",")))
    parms$wbar = as.numeric(unlist(strsplit(input$wbar,",")))
    parms$cbar = as.numeric(unlist(strsplit(input$cbar,",")))
    parms$alpha = as.numeric(unlist(strsplit(input$alpha,",")))
    parms$eta = as.numeric(unlist(strsplit(input$eta,",")))
    parms$p = as.numeric(unlist(strsplit(input$p,",")))
    parms$q = as.numeric(unlist(strsplit(input$q,",")))
    parms$s = as.numeric(unlist(strsplit(input$s,",")))
    parms$w = input$w
    parms$r = input$r
    parms$m = as.numeric(unlist(strsplit(input$m,",")))
    parms$it = input$it
      
    Y = do.call(simulatey,parms)
    
    goods = c("Corn","Iron","Sugar")
    
    olf = Y$agg_employ*Y$w
    bolf = rowSums(cbind(Y$profits1,Y$profits2,Y$profits3))
    golf = Y$agg_demand*(Y$r/(1+Y$r))
    dolf = rowSums(cbind(olf,bolf,golf))
    
    
    suppressWarnings(grid.arrange(
                 graffy(cbind(Y$p1,Y$p2,Y$p3),goods,title="Prices",leg=TRUE),
                 graffyrex(Y),
                 graffy(cbind(Y$q1,Y$q2,Y$q3),goods,title="Quantities"),
                 graffy(cbind(Y$s1,Y$s2,Y$s3),goods,title="Inventories"),
                 graffy(cbind(Y$costs1,Y$costs2,Y$costs3),goods,title="Sectoral costs"), 
                 graffy(cbind(Y$profits1,Y$profits2,Y$profits3),goods,title="Sectoral profits"),
                 graffy(cbind(Y$employments1,Y$employments2,Y$employments3),goods,title="Sectoral employment"),
                 graffy(cbind(Y$h1,Y$h2,Y$h3),goods,title="Costs per unit produced"),
                 graffy(cbind(olf/dolf,bolf/dolf,golf/dolf),c("Wages","Profits","Interest"),title="Income distribution share",leg=TRUE),
                 graffy(cbind(Y$m1,Y$m2),c("Workers","Capitalists"),title="Savings",leg=TRUE),
                 dumbgraffy(bolf,"Total profit"),
                 dumbgraffy(Y$agg_employ,"Aggregate employment"),
                 dumbgraffy(Y$agg_demand,"Aggregate demand"),
                 dumbgraffy(Y$w,"Wage rate"),
                 dumbgraffy(Y$r,"Interest rate"),
                 nrow = 5))
    
  })
}

ui <- fluidPage(
  titlePanel("Nonlinear Dynamic Model of Classical Macrodynamics from Chapter 7 of Ian Wright's Ph.D. Thesis", windowTitle = "Hello Shiny!"),
  titlePanel(tags$h4(
    tags$a("The Law of Value: A Contribution to the Classical Approach to Economic Analysis",href="http://pinguet.free.fr/wrightthesis.pdf")
  )
  ),
        tabsetPanel(
          tabPanel("Model specification",
             fluidRow(
               column(3,
      textInput("A1",
                label = "Corn input per unit (Corn, Iron, Sugar)",
                value = paste(as.character(parms$A[1,]), collapse=",")
                ),
      textInput("A2",
                  label = "Iron input per unit (Corn, Iron, Sugar)",
                  value = paste(as.character(parms$A[2,]), collapse=",")
                ),
      textInput("A3",
                  label = "Sugar input per unit (Corn, Iron, Sugar)",
                  value = paste(as.character(parms$A[3,]), collapse=",")
                ),
      textInput("l",
                label = "Labor input per unit (Corn, Iron, Sugar)",
                value = paste(as.character(parms$l), collapse=",")
                ),
      textInput("eta",
                label = "Elasticities (3 Prices*, 3 Quantities*, Wages, Interest)",
                value = paste(as.character(parms$eta), collapse=",")
                )),
               column(3,
      textInput("wbar",
                label = "Worker consumption bundle (Corn, Iron, Sugar)",
                value = paste(as.character(parms$wbar), collapse=",")
                ),      
      textInput("cbar",
                label = "Capitalist consumption bundle (Corn, Iron, Sugar)",
                value = paste(as.character(parms$cbar), collapse=",")
                ),
      textInput("alpha",
                label = "Worker and capitalist propensity to spend",
                value = paste(as.character(parms$alpha), collapse=",")
                ),
      textInput("m",
                label = "Money distribution at t0 (Workers,Capitalists)*",
                value = paste(as.character(parms$m), collapse=",")
                )),
            column(3,
      textInput("p",
                label = "Prices at t0 (Corn, Iron, Sugar)*",
                value = paste(as.character(parms$p), collapse=",")
                ),
      textInput("q",
                label = "Quantities at t0 (Corn, Iron, Sugar)",
                value = paste(as.character(parms$q), collapse=",")
                ),
      textInput("s",
                label = "Inventories at t0 (Corn, Iron, Sugar)*",
                value = paste(as.character(parms$s), collapse=",")
                ),
      numericInput("w",
                label = "Wage rate at t0",
                value = parms$w
                ),
      numericInput("r",
                label = "Interest rate at t0",
                value = parms$r
                )),
                column(3,
      numericInput("it",
                label = "Iterations*",
                value = parms$it
                ))
    
    )  
    ), tabPanel("Results",
              fluidRow(
                column(12,
                    plotOutput("everything",height = 1300)
                      )
                      )
              )
  )
)

shinyApp(ui = ui, server = server)