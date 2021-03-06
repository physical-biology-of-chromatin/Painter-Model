# Painter-Model
A unified quantitative framework to systematically characterize epigenome regulation and memory 
Uses Gillespie algorithm to simulate. Ref https://pubs.acs.org/doi/10.1021/j100540a008

Run the Painter_Model.py.  
## Parameters
  >_K0_                  Nucleosome turnover rate\
  >_Kme_                 Acetylation, Methylation rate\
  >_Eme_                 Epsilon - spreading efficiency\
  >_r_                   Allosteric boost\
  >_delta_               Reader Writer Mechanism\
  >_gamma_               Contact Probability |i - j| exponent\
  >_delT_                Determines the time interval between recording epigenetic configuration state<br/> 
  >_limt_                Number of Trajectories or cells<br/> 
  >_ts_                  Specifies the time at which painter unbinds  
  >_a0_                  Maximal transcription rate\
  >_d_                   Steepness of transcriptional switch<br/> 
  >_mc_                  Critical P(M) value at which the switch happens<br/>
  >_m_decay_             mRNA decay rate\
  >_g_decay_             protein decay rate\
  >_m_gfp_               mRNA translation rate\
  >_g_mat_               protein maturation rate\
  >_{i_start,i_end}_     specifying the painter region
  
For transcription, transcribing region can be specified in line 162<br/>
For enzyme limitation uncomment line 132. 
  >_Ntot_        Total number of HMEs\
  >_rconst_      ku/kb = ratio binding-unbinding rates

(Use appropriate parameter values mentioned in the article to obtain the results)  
 
## Output
   Output is stored in the variable fin_out. Line 226 specifies transcription or chromatin state output to be stored. 
   'Output.txt' stores single trajectory data as well.  
   
