#### Server ####
server <- function(input, output) {
  # Source: https://github.com/open-meta/uiStub
  cat("Session started.\n") # this prints when a session starts
  onSessionEnded(function() {cat("Session ended.\n\n")})  # this prints when a session ends
  
  
  ############# PAGE 1 ##################
  #Define Variables
  R <- reactive({input$R})
  S <- reactive({input$S})
  h2 <- reactive({R()/S()})
  popsize <- reactive({input$popsize})
  
  ##################Create Dataset#################
  # Assume traits are 
  gen1 <- reactive({
    rnorm(popsize(), 0, 1)
  })
  gen2 <- reactive({
    rnorm(popsize(),R(),1)
  })
  dat <- reactive({
    data.frame(tot = c(gen1(),gen2()),
               gen=c(rep("Gen1",popsize()),rep("Gen2",popsize())))
  })
  
  mu <- reactive({
    ddply(dat(), "gen", summarise, grp.mean=mean(tot))
  })
  
  output$distPlot <- renderPlot({
    # Validation source: 
    # https://shiny.rstudio.com/reference/shiny/0.14/validate.html
    validate(
      need(input$R <= input$S, "Usage Error: R > S")
    )
    # Hist source: 
    # http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization? 
    ggplot(dat(), aes(x=tot, color = gen)) + 
      geom_histogram(fill="white", binwidth=.1) +
      facet_grid(gen ~ .) + 
      labs(title="Population Distribution over 2 Generations",
           subtitle=paste("N = ", popsize(),"; h^2 = ", h2()),
           x="Trait", y= "Frequency") +
      scale_color_manual(values=c("purple","#1da3a5")) +
      geom_vline(data=mu(), aes(xintercept=grp.mean, color=gen),
                 linetype="dashed") +
      theme_bw() 
    
  })
  

  ############# PAGE 2 ##################
  #Define Variables
  threshhold <- reactive({input$threshhold})
  N <- reactive({input$N})
  deltaN <- reactive({input$deltaN})
  pop <- reactive({randpeople(N())}) 
  s <- reactive({input$strength})
  # NS acts on pop
  selpop <- reactive({naturalselection(pop(),threshhold(),s())})
  # recombine pop
  newpop <- reactive({recombinepop(selpop(),deltaN())})

  ##################Create Dataset#################
  output$gen1Plot <- renderPlot({
    validate(
      need(N() > 0, "All was quiet on the broken simulation.")
    )
    plotheight(pop())
  })
  
  output$gen2Plot <- renderPlot({
    validate(
      need(deltaN() > 0, "All was quiet on the broken simulation.")
    )
    validate(
      need(threshhold() < 8, "Everyone died! This isn't a very interesting simulation at all anymore...")
    )
    plotheight(newpop())
  })

  output$tblgen1 <- renderDataTable({
    validate(
      need(N() > 0, "All was quiet on the broken simulation.")
    )
    countfreq(pop())
  })
  
  output$tblgen2 <- renderDataTable({
    validate(
      need(deltaN() > 0, "All was quiet on the broken simulation.")
    )
    validate(
      need(threshhold() < 8, "Congrats, you killed everyone... I hope you're happy.")
    )
    countfreq(newpop())
  })

  
  
  ########################### PAGE 3 ###############################
  # Define variables
  NTrials <- reactive({input$NTrials})
  Ne <- reactive({input$Ne})
  Ngen <- reactive({input$Ngen})
  fAA <- reactive({input$fAA}) 
  fAa <- reactive({input$fAa}) 
  
  # Run simulation
  sims <- reactive({
    dat <- fixitysimulator(NTrials(),Ne(),Ngen(),fAA(),fAa())
    rbind(dat, c("Means", mean(dat$ReachFix),round(mean(dat$GensToFixity, na.rm =T), digits=2)))
  })
  
  output$tblneutral <- renderDataTable({
    datatable(sims())
  })
  
  output$textvar <- renderText({ 
    paste("Allele frequency (freq(A) or p) =", sprintf("%.2f",round(fAA() + .5 * fAa(), digits=2)))
  })
  
  output$aggplot <- renderPlotly({
    aggplot
  })
  
  output$allele10pop <- renderPlotly({
    allele10pop
  })
  
  output$allele10time <- renderPlotly({
    allele10time
  })

}
