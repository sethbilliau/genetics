library(shiny)
library(shinythemes)
library(ggplot2)
library(plyr)
library(shinyLP)
library(DT)
library(plotly)
source("GENETICS.R")
source("NEUTRAL.R")
source("KIMURA.R")

#tags: https://shiny.rstudio.com/articles/tag-glossary.html

cat("Application started...\n")
#### UI ####
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
     .sel {
        margin-left:25px;
        list-style-type: disc;
        list-style-position: inside;
        text-indent: -1em;
     }
      
      img {
        max-width:100%;
        max-height:100%;
      }

                    "))
    ),
  # Source: https://rstudio.github.io/shinythemes/ 
  theme = shinytheme("cerulean"),
  
  
  
  ############## HOME PAGE ####################
  navbarPage("Simulating Population Genetics",
    tabPanel("Home",
      fluidRow(
        column(12, align="center",
               img(src="darwin.jpg", height=400)
               # Source: National Portrait Gallery, London, courtesy of A. Berry 
               # https://www.summitpost.org/r-a-fisher/1022093
        )
      ),
      h1("Simulating Population Genetics with RShiny", align="center"),
      # Source: https://loremipsum.io/
      p("My final project explores population genetics through a series of simulations. In the first section, we will introduce a simple and interactive model of directional selection, later improving upon our model by introducing Fisherian polygenetic inheritance with a similar interactive simulation.",
        "The second section explores how long it takes for a neutral mutation to reach fixity with a given allele frequency and population size, following in the footsteps of Kimura and Ohta.",
        "The goal of this project is to dig deeper into results from class, using simulation to illustrate interesting results, answer complex questions, and verify analytical solutions.", align="center"),
      hr(),
      p("Image Source: C. Darwin, Wikipedia Creative Commons,", tags$a(href="https://commons.wikimedia.org/wiki/File:Charles_Robert_Darwin_by_John_Collier.jpg", "URL link."))
    ),
    
    
    
    ############## PAGE 1 ####################  
    navbarMenu("Selection", 
      tabPanel("Simple Directional Selection",
        # App title ----
        titlePanel("A Simple Illustration of Directional Selection"),
        fluidRow(
          column(8,
                 h4("Purpose of the Simulation:"),
                 p("In", tags$em("The Origin of Species,"), "Charles Darwin articulated his theory of evolution by natural selection, hypothesizing that favored, heritable traits are likely to be",
                   "preserved in future generations while unfavorable traits are likely to be eliminated. One of the more interesting forms of natural selection is directional selection. Formalized by August Weismann", 
                   " and Ivan Schmalhausen,",
                   "directional selection occurs when an extreme phenotype is favored, leading to a shift in the mean phenotypic distribution of a population over time. We will begin with a very basic simulation describing directional selection."),
                 h4("Explaining the Simulation:"),
                 p("This simulation prompts the user for three inputs: population size, selection differential and response to selection. Definitions for each input are provided.",                   
                   "If you play around a bit with the simulation, it becomes clear that as the response to selection",
                   "increases, the distribution of the second generation moves to the right - directional selection at work. We can also see that ",
                   "as the selection differential increases, heritability of a given trait decreases given that R is constant. This is a good, basic model for directional selection and a good reminder of the relationship between heritability, response to selection, and the selection coefficient (h^2 =R / S).")
          ),
          column(4,
                 img(src="schmal.jpeg", height=400),
                 p("Ivan Schmalhausen")
          )
        ),
        tags$hr(),
          # Sidebar layout with input and output definitions ----
          sidebarLayout(
           # Sidebar panel for inputs ----
           sidebarPanel(
             # Input: Numerics for the popsize, S and R ----
             # Source: https://shiny.rstudio.com/reference/shiny/latest/numericInput.html
             numericInput(inputId = "popsize",
                          label = "Population Size", 
                          value = 10000, min = 1, max = 10^6, step = 100,
                          width = NULL),
             p("Larger Population sizes will reduce the variability of the simulation"),
             numericInput(inputId = "S",
                          label = "Selection Differential (S):", 
                          value = 1, min = 0, max = 6, step = .01,
                          width = NULL),
             p("The Selection Differential is the difference in means between the Gen. 1 total population and the subset of the Gen. 1 population favored by selection."),
             numericInput(inputId = "R",
                          label = "Response to Selection (R):", 
                          value = .5, min = 0, max = 6, step = .01,
                          width = NULL),
             p("The Response to Selection is the difference in means between the Gen. 1 total population and the subset of the Gen. 2 total population.")
           ),
           
           # Main panel for displaying outputs ----
           mainPanel(
             
             # Output: ggplot Histogram ----
             plotOutput(outputId = "distPlot")
           )
         ),
         # Footer Source: 
         # https://stackoverflow.com/questions/30205034/shiny-layout-how-to-add-footer-disclaimer
         hr(),
          h4("Criticisms:"),
          p("This model makes a number of assumptions to simplify our simulation. Here is a short (nonexhaustive) list of some things that should give us pause:"),
          tags$li("Both populations are generated from a normal distribution with an variance of 1. The means differ based on R, but the variances do not. It is entirely possible that selection could alter the variance of a trait between populations, so this assumption of equal variances may not be a good one.", class ="sel"),
          tags$li("Because this model generates both generation populations from separate distributions, there's not actually any simulation of heredity or transfer of genetic material. Though this does illustrate how these processes would work, it doesn't actually use alleles from previous generations to form new generations.", class = "sel"),
          tags$br(),
          p("Our next simulation will attempt to address these concerns. Navigate to 'Directional Selection with Polygenics' to continue."),
         hr(),
         p("Inspired by Dr. Berry's Lecture 10 Notes, see my bibiography for all content-related citations.", tags$a(href= "https://link.springer.com/article/10.1134%2FS1019331609030149", "Image Credit."))
           
      ),
      
      
      
      
      ############################ PAGE 2 ##################################
      tabPanel("Directional Selection with Polygenics",
       titlePanel("Illustrating Natural Selection: Height determined by Four Loci"),
       fluidRow(
         column(8,
            h4("Purpose of the simulation:"),
            p("In his ground-breaking paper titled", tags$a(href="http://www.indiana.edu/~curtweb/L567/readings/Fisher1918.pdf", "'The Correlation between Relatives on the Supposition of Mendelian Inheritance',"), "geneticist and statistician R.A. Fisher formalized the first theory of polygenetic inheritance, reconciling Mendelian (discrete) and biometrical (continuous) genetics.", 
              "His basic idea was that many genetic loci make small, discrete contributions to determine a single phenotypic characteristic.",
              "Though each gene is discrete, this results in a continuous, normally distributed phenotypic range. Fisher also notes that only a small number of genes are required to yield a normally distributed phenotypic distribution ('The Correlation between Relatives on the Supposition of Mendelian Inheritance').",
              "This simulation will refine our previous simulation describing directional selection, using polygenetic inheritance as our method of producing continuous phenotypic variation."),
            h4("Explaining the simulation:"),
            p("This simulation works by generating a user-specified number of Generation 1 individuals, each with 4 genetic loci. Like our example in lecture, each gene has 2 alleles. One is capital letter and the other is a lowercase letter yielding 3 unique genotypes (ex. 1 locus has the genotypes AA, Aa, aa). Each individual is generated with each genotype equally likely. Each capital letter contributes one unit height, yielding 9 distinct height phenotypes.",
              " The phenotypic and genotypic frequency distributions for Generation 1 are visualized for the user.",
              " As R.A. Fisher's work would predict, the distribution of phenotypes is approximately normally distributed."),
            p("Here's where natural selection kicks in. In this scenario, individuals are deemed unfit if they fall beneath the user specified height threshhold. ",
              "The selection coefficient, the percentage of unfit individuals that are unable to reproduce, is specified by the user.", 
              " Then, a user-specified number of Generation 2 individuals are generated from the set of individuals able to reproduce.",
              " These individuals are generated in a traditional Mendelian manner, each parent passing on one of their allele's to a generation 2 child. Similar plots and tables are displayed below those for generation 1.")
          ),
         column(4,
                img(src="fisher.jpg", height=400),
                p("R.A. Fisher")
         )
        ),
        tags$hr(),
       # Sidebar layout with input and output definitions ----
       sidebarLayout(
         
         
         # Sidebar panel for inputs ----
         sidebarPanel(
           # Input: Numerics for the popsize, S and R ----
           # Source: https://shiny.rstudio.com/reference/shiny/latest/numericInput.html
           p("Choose population sizes for generations 1 & 2 below:"),
           numericInput(inputId = "N",
                        label = "Gen 1 Population Size", 
                        value = 1000, min = 1, max = 10^6, step = 100,
                        width = NULL),
           numericInput(inputId = "deltaN",
                        label = "Gen 2 Population Size", 
                        value = 1000, min = 1, max = 10^6, step = 100,
                        width = NULL),
           sliderInput(inputId = "threshhold",
                       label = "Height threshhold for Survival:", 
                       value = 1, min = 0, max = 8, step = 1,
                       width = NULL),
           p("Note: In our simulation, natural selection acts against all individuals less than the specified height threshhold. 'Fit' individuals are taller than the threshhold and 'unfit' individuals are equal to or below it. (ex. a height threshhold of 1 means that individuals of one height unit or less are 'unfit')"),
           sliderInput(inputId = "strength",
                       label = "Selection Coefficient (S):", 
                       value = .5, min = 0, max = 1, step = .01,
                       width = NULL),
           p("Note: This is the percentage of unfit individuals that are unable to pass on their genes due to natural selection. (ex. S = 1 prohibits all unfit individuals from passing on their genes)")
        ),
  
         
         # Main panel for displaying outputs ----
         mainPanel(
           # Output: ggplot Histogram ----
           h2("Generation 1 Phenotypic Distribution"),
           plotOutput(outputId = "gen1Plot"),
           h2("Generation 1 Genotypic Distribution"),
           dataTableOutput('tblgen1'),
           p("Note: By counting the number of capital letters, we establish genotype frequencies."),
           h2("Generation 2 Phenotypic Distribution"),
           plotOutput(outputId = "gen2Plot"),
           h2("Generation 2 Genotypic Distribution"),
           dataTableOutput('tblgen2'),
           p("Note: Tables will likely be inaccurate if population sizes are small (< 10). It's a bug I'm working on.")
         )
       ),
       hr(),
       h4("Criticisms:"),
       p("Clearly, this simulation is much more elegant than our first attempt. However, we could add more complexity to our simulation by adding:"),
       tags$li("Sexual Selection working opposite Natural Selection", class ="sel"),
       tags$li("Dominant and Recessive alleles", class = "sel"),
       tags$li("Relationships between loci", class = "sel"),
       tags$li("Environmental Factors", class = "sel"),
       tags$br(),
       p("Given more time, I will add these to another page in the simulation. For now, I will move on from simply illustrating results from class and begin trying to answer genetic questions difficult to answer analytically using methods of simulation."),
       # Footer Source: 
       # https://stackoverflow.com/questions/30205034/shiny-layout-how-to-add-footer-disclaimer
       hr(),
       p("Inspired by Dr. Berry's Lecture 11 Notes, see my bibiography for all content-related citations,", tags$a(href= "https://link.springer.com/article/10.1134%2FS1019331609030149", "Image Credit."))
               
      )
    ),
    
    navbarMenu("Neutral Theory",
      tabPanel("Neutral Fixation - Interactive Simulation",
         titlePanel("How long does it take for a neutral allele to reach Fixity?"),
         fluidRow(
           column(8,
                  h4("Purpose of the Simulation: "),
                  p("The purpose of this simulation is to test how long it takes a neutral allele A to reach fixity. The mechanics of this simulation are taken from knowledge of neutral alleles", 
                            "and fixation given in lecture and from a reading of Motoo Kimura's 1962 paper entitled 'On the probability of fixation of mutant genes in a population'."),
                  h4("Explaining the simulation:"),
                  tags$li(p(tags$b("Set-up:"), "The number of generations required to reach to fixity will depend on the initial allele frequency of A and the effective population size of our species. ",
                    "These constraints will be specified by the user, and a generation 0 of individuals will be created according to these initial conditions."), class ="sel"),
                  tags$li(p(tags$b("Describing a Trial:"), "The simulation will simulate random mating between individuals, keeping effective population size constant while permitting genetic recombination.",
                    "After mating, the simulation will check to see if our allele A has reached 100% frequency. If this condition is met, the simulation will record how many generations were required.",
                    "If the condition is not met, the process will repeat to with the new generation. This procedure will",
                    " continue until allele A reaches fixity and the number of required generations is recorded or the user-specified maximum number of generations is reached. This process describes one trial."), class ="sel"),
                  tags$li(p(tags$b("Data Visualization:"), "The user-specified number of trials are performed. The results of each trial are recorded in a datatable. Mean statistics are given in the final row of the data table.", 
                                 "The Mean of the 'ReachFix' category is the proportion of trials that reached fixity in the given time frame. We would expect this to be approximately equal to the initial allele frequency of our neutral allele A.", 
                                  "The mean of 'GensToFixity' is the arithmetic mean of generations that actually reach fixation ignoring generations in which the allele is lost. You could argue this is undesirable (ie. we ignore large values of GensToFixity since", 
                                 "we terminate the simulation after a specified number of generations.).",
                                 "From a population genetics perspective, however, this value is of great interest when studying gene substitution as we will see in the second section of 'Neutral Theory'."), class ="sel")
                  ),
           column(4,
                  img(src="kimura.jpg", height=400),
                  p("Motoo Kimura")
           )
         ),
         tags$hr(),
              # Sidebar layout with input and output definitions ----
               sidebarLayout(
                 
                 
                 # Sidebar panel for inputs ----
                 sidebarPanel(
                   # Input: Numerics for the popsize, S and R ----
                   # Source: https://shiny.rstudio.com/reference/shiny/latest/numericInput.html
                   h4("Biological Parameters:"),
                   numericInput(inputId = "Ne",
                                label = "Effective Population Size (Ne):", 
                                value = 15, min = 1, max = 100,
                                width = NULL),
                   p("Note: Recall that the effective population size is defined as the number of
                     individuals in a theoretically ideal population having the
                     same magnitude of genetic drift as the actual population."),
                   sliderInput(inputId = "fAA",
                               label = "Frequency of homozygous dominant (freq(AA)):", 
                               value = .5, min = 0, max = 1, step = .01,
                               width = NULL),
                   sliderInput(inputId = "fAa",
                               label = "Frequency of heterozygote (freq(Aa)):", 
                               value = .5, min = 0, max = 1, step = .01,
                               width = NULL),
                   tags$b(textOutput("textvar")),
                   tags$hr(),
                   h4("Simulation Constraints:"),
                   p("These constraints are artificial, but they prevent stack overflow and runtime"),
                   numericInput(inputId = "NTrials",
                                label = "Number of Simulations:", 
                                value = 10, min = 1, max = 100, step = 1,
                                width = NULL),
                   p("Note: Increasing this value will give you more data points, but it will take longer to run."),
                   numericInput(inputId = "Ngen",
                                label = "Maximum number of Generations per Simulation:", 
                                value = 100, min = 0, max = 100, step = 1,
                                width = NULL)
                   ),
                 
                 
                 # Main panel for displaying outputs ----
                 mainPanel(
                   # TODO
                   dataTableOutput('tblneutral')
                  )
               ),
               hr(),
               h4("Criticisms:"),
               p("All in all, I'm really proud of this simulation!", 
                 "Here are a few one way I think it could be improved: "),
               tags$li(p(tags$b("Time Complexity and RAM trouble:"), 
                         "I have not optimized my code in regards to memory allocation.",
                         "This could result in stack overflow for large input values. Even if RAM is not the issue, I have not optimized my algorithm for time complexity either. 
                          Therefore, large input values may take a long time to run.", 
                         "I think this interactive UI app is valuable as it gives the user an interactive way to explore and construct simple examples of a neutral allele going to fixity.",
                         "However, it's not efficient to use this version of the app for analysis. For these reasons, I'm",
                         "going to work with a modified version of the app's backend, nixing the responsive UI for the sake of time and RAM in the next section."), class ="sel"),
               p("Please continue to the next tab in the 'Neutral Theory' section of the navbar."),
               # Footer Source: 
               # https://stackoverflow.com/questions/30205034/shiny-layout-how-to-add-footer-disclaimer
               hr(),
               p("Inspired by Dr. Berry's Lecture 13-14 Notes, see my bibiography for all content-related citations,", tags$a(href="http://evolutionlist.blogspot.com/2006/03/", "Image Credit."))

     ),
     tabPanel("Testing Kimura-Ohta: Generations to Fixity",
              titlePanel(paste("Kimura-Ohta on generations to fixity")),
              fluidRow(
                column(8,
                       h4("Purpose of the simulation:"),
                       p("In a paper entitled", a(href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1212239/pdf/763.pdf", "'The Average Number of Generations until Fixation of a Mutant Gene in a Finite Population'"),
                                 "written by Motoo Kimura and Tomoko Ohta, the authors demonstrate that as the initial allele frequency of a given neutral allele A goes to zero,",
                                 "the number of generations it takes for A to reach fixation is four times the effective population size. The analytical derivations in the paper are extremely elegant. 
                                  However, they performed all of their analysis on two very old computers, the TOSBAC 3400 (this",tags$a(href="http://museum.ipsj.or.jp/en/computer/main/0004.html", "computer's speed"), "is unknown)", 
                                  " and the IBM 360 (the",  
                                  tags$a(href="https://www.ibm.com/ibm/history/exhibits/mainframe/mainframe_PP2091.html", "1967 model"), "could likely perform 16.6 MIPS or million instructions per second).",
                                  "Because of these constraints, the scientists only performed 400 trials for each permutation of 2 effective population sizes (Ne = 5, 10) and 9 starting allele frequencies (.1 to .9 by steps of .1 for Ne = 5, 0.05 to 0.95 by steps of .05 for Ne = 10).", 
                                 "In this simulation, I hope to use modern computers to verify the Kimura-Ohta result, increasing the trial sizes to decrease sampling error and to test a wider variety of population sizes.", 
                                 "I will also record how long each function takes to execute on a modern computer (1.2 GHz processor - approximately 1200 MIPS) in an attempt to compare the difficulty of performing this simulation in two eras.")
                ),
                column(4, 
                       img(src="ohta.png"),
                       p("Tomoka Ohta")
                )
              ),
              tags$hr(),
              h2("Running the Simulation:"),
              p("After making a few alterations to the backend of the app on the previous page, I ran simulations an initial allele frequency", 
                "of 0.1 and population sizes varying from 10 through 100 by steps of 10. In each population size, I generated 1000 trials. The figure below gives my results. Each point represents 1000 trials with one effective population size.",
                "The red line gives the Kimura-Ohta prediction for generations to fixity as p goes to 0. The blue line gives a linear regression of my observed points with a shaded region for standard error. Note that they are reasonably close, even though p = 0.1 is not all that",
                "close to 0. Note that all of these graphs are interactive! Feel free to play around with them, zooming in different places to make things clear."),
              plotlyOutput("allele10pop"),
              p("To my surprise, this simulation took over 5 hours to run. After studying my code and comparing a few regression coefficients, I hypothesized that",
                "the time it takes for my code to complete likely increases quadratically as effective population size increases linearly. This is hypothesis is supported by the R^2 values of exponential (0.3745), linear (0.9414),",
                "and quadratic (0.9938) regression models for my observed data. I have plotted time against Ne in the graph below with a red line representing quadratic regression. Extrapolating point estimates from the quadratic model, it would take around 2.65, 7.07, and 14.5 hours to do the same simulation with",
                "effective population sizes of 200, 500, and 1000 respectively. Even my modern computer with a processor around 72 times faster than the IBM 360, this simulation takes a long time to complete!",
                "This certainly gave me a greater appreciation for the Kimura-Ohta paper."),
              plotlyOutput("allele10time"),
              p("I continued my experiment by testing all 45 permutations of nine different allele frequencies (from 0.1 to 0.9 by steps of 0.1) and 5 population sizes (from 10 to 50 by steps of 10). The results are plotted below. I also plotted the predicted analytical results given in equation 12 of the Kimura-Ohta paper",
                "a corresponding color. Note though the correspondence is somewhat satisfactory though it seems the Kimura-Ohta analytical result almost always overestimates the number of generations to reach fixity. This result also occurred in the original paper though it is largely unexplained. This graphic is a bit difficult to read",
                "as originally rendered, so please make use of the graph's ability to scale and zoom."),
              plotlyOutput("aggplot"),
              p("All in all, the Kimura-Ohta paper is a remarkable scientific achievement, especially in the context of its era. I'm really rather proud of these graphics and this result as I think it does a good job of expanding upon the work of ",
                "past scientists."),
              hr(),
              p("Inspired by Dr. Berry's Lecture 13-14 Notes, see my bibiography for all content-related citations,", tags$a(href="https://en.wikipedia.org/wiki/Tomoko_Ohta#/media/File:Tomoko_Harada_cropped_1_Tomoko_Harada_201611.png", "image Credit."))
              
    )
              
    ),
  tabPanel("Bibliography/Video", 
           tags$a(href = "https://docs.google.com/document/d/1p1I-np_XhgQBpXMOoOCrl6NiHe0lB4w2JuiCLqTnE_A/edit?usp=sharing", "Click Here to view my Bibliography for all things Population Genetics."),
           tags$br(),
           tags$a(href = "https://www.youtube.com/watch?v=5DExmABnXwc&feature=youtu.be", "Click Here to watch a YouTube video of me explaining what the app does (it's super long and boring).")
  )
  ),
  
  cat("UI Successfully Loaded...\n")
)
