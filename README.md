# Painter-Model
A unified quantitative framework to systematically characterize epigenome regulation and memory 
Uses Gillespie algorithm to simulate. Ref https://pubs.acs.org/doi/10.1021/j100540a008

Give in the appropriate parameter values mentioned in the article. 
## Parameters
  K0                  Nucleosome turnover rate
  Kme                 Acetylation, Methylation rate
  Eme                 Epsilon - spreading efficiency
  r                   Allosteric boost
  delta               Reader Writer Mechanism
  gamma               Contact Probability |i - j| exponent
  delT                Determines the time interval between recording epigenetic configuration state. (Time resolution) 
  limt                Number of Trajectories or cells 
  ts                  Specifies the time at which painter unbinds    
  #TRANSCRIPTION
  a0                  Maximal transcription rate. Actual would be higher considering leaky rate 
  d                   Steepness of transcriptional switch  
  mc                  Critical P(M) value at which the switch happens
  m_decay             mRNA decay rate
  g_decay             protein decay rate
  m_gfp               mRNA translation rate
  g_mat               protein maturation rate
  {i_start,i_end}     specifying the painter region
  
For transcription transcribing region can be specified in line 122. 

## Output
   Output is stored in the variable fin_out. In line 186 transcription or chromatin state results can be specified as needed. 
   'Output.txt' stores a single tajectory data as well.  
