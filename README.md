# 2023 SPWLA PDDA Machine Learning Competition
## Automatic Well-Log Depth Shift with Data-Driven Methods

<img style="float: left;" src="image/pdda.png" alt="pdda" title="" width="300" height="100"/>
<!-- 
- <a>Call for Registration </a>
- <a>Sponsoring Opportunities</a>
- <a>Challenge Description</a>
- <a>Evaluation</a>
- <a>Timeline</a>
- <a>Competition Rules</a>
- <a>Contest Committee</a>
-->

### <a>Call for Registration</a>
Dear petrophysics and data science enthusiasts,

We are thrilled to announce our upcoming **Automatic Well-Log Depth Shift with Data-Driven Methods** machine learning competition, designed to challenge participants to develop innovative solutions for the critical task of depth-shift correction in well logs.

Well logs are crucial data sources for the oil and gas industry, providing critical information about subsurface formations and geological properties. However, the data collected through well logs can often be affected by borehole environments, such as the sticking and slipping of logging tools, resulting in inaccurate and unreliable results. Manual depth shift of well logs is time-consuming and subject to the expertise of interpreters. To automate this process and improve the accuracy of depth shift, this is where your skills in machine learning can make a significant difference.

We invite all petrophysicists, data scientists, machine learning engineers, and anyone passionate about data-driven solutions to join our competition. The objective is to develop a depth shift correction model that can accurately and efficiently correct depth shifts in well logs.

The competition will be hosted on a Codalab and will run for eight weeks. Participants will have access to a well-curated dataset, including a training set of 9 depth-aligned well log set and a test set of 3 misasligned well log set. They will develop data-driven model to automatically predict the depth shift for misaligned well logs in the test set. Their solutions will be evaluated based on the prediction accuracy of depth shift and the shifted well-log values. More details can be found in the tutorial provided.

The winners will receive attractive prizes, including non-cash rewards(could be gift cards), certificates, and a chance to showcase their solution to industry experts in the upcoming SPWLA conferences.

So, if you are interested in testing your machine learning skills and contributing to the oil and gas industry, register now and join us in this exciting competition.

To register, please send an email to **pdda_sig@spwla.org** with your team name, team members, and the email address of the team lead.

We look forward to seeing your innovative solutions.

Best regards,

Wen Pan


### <a>Sponsoring Opportunities</a>
Attention all companies interested in supporting data-driven solutions and innovative technologies!

We are excited to announce our upcoming **Automatic Well-Log Depth Shift with Data-Driven Methods** machine learning competition and are seeking sponsors to help make it a success.

This competition aims to challenge participants to develop innovative solutions for automatic depth shift correction for well logs, which are crucial data sources for the oil and gas industry. By sponsoring this competition, your company can help us promote this critical task and encourage data science enthusiasts to develop solutions that will help advance the industry.

Sponsoring our competition is an excellent opportunity for your company to showcase your brand to a targeted audience of petrophysicists, petroleum engineers, data scientists, machine learning engineers, and other professionals interested in data-driven solutions. You will receive prominent visibility on our website, conferences and social media channels, and we will acknowledge your support in our press releases (journals and conference preceedings) and other promotional materials.

We invite all interested companies to consider sponsoring our event. Please contact us via our email **pdda_sig@spwla.org**.

Thank you for considering our proposal. We look forward to hearing from you soon.

Best regards,

Wen Pan

### <a>Challenge Description</a>
The objective of this contest is to develop a data-driven model to shift misaligned well logs to a reference GR log. Well logs of nine wells from the same field that have been depth shifted by a petrophysicist. Participants' models will be trained/validated with these nine wells and tested on misaligned well logs from other three wells from the same field.

The submission should include the shifted well logs and the shifted depths



### <a>Evaluation and Scoring Website</a>
Please register your team at Codalab with your team leader's email. A link to the competition will be sent to you to access the scoring system.
You will be given nine wells with well logs aligned, then you will train a model to perform depth shift for three hidden wells with misaligned well logs. The gamma ray log (GR) is used as the reference log, you are required to shift the other well logs (Density (RHOB), Neutron (NPHI), and Resistivity (RD)) to align with the reference GR log. You will predict the values of corrected well logs (RHOB_pred, NPHI_pred, RD_pred), and the asssocited depth shift (RHOB_dept_pred, NPHI_dept_pred, RD_dept_pred) for the three test wells. You will be scored based on the normalized mean squared error (NMSE) of your well-log prediction and the Mean absolute deviation (MAD) of your depth shift prediction.

Submissions are evaluated based on normalized mean squared error(NMSE) of corrected well logs and mean absolute deviation (MAD) of depth shift prediction, the final score will be rank transformed and averaged to avoid the scaling of different metrics.

The metrics for scoring are as follows:

$$NMSE = \frac{1}{mn}\frac{1}{Var(\mathbf{y})}\sum_{i=1}^{m}\sum_{j=1}^{n}(\hat{\mathbf{y_{i,j}}} - \mathbf{y_{i,j}})^{2}$$

where
- $\hat{y_{i,j}}$ the prediction (RHOB_pred, NPHI_pred, RD_pred) of the **values** of shifted well log j for sample i, $y_{i,j}$ is the actual **values** of the well log j for sample i shifted by a petrophysicist. 
- $m$ is the sample size.
- $n$ is the number of well logs (RHOB, NPHI, RD): 3.
- $Var$ is the variance.

$$MAD = \frac{1}{mn}\sum_{i=1}^{m}\sum_{j=1}^{n}|\hat{\mathbf{d_{i,j}}} - \mathbf{d_{i,j}}| $$

where
- $\hat{d_{i,j}}$ is the prediction (RHOB_dept_pred, NPHI_dept_pred, RD_dept_pred) of **depth shift** of raw well log j for sample i, $d_{i,j}$ is the actual **depth shift of raw well logs** by a petrophysicist. 
- $m$ is the sample size.
- $n$ is the number of well logs (RHOB, NPHI, RD): 3.


**Note**:
- Please remember to use random_state for all randomization steps, so the results are reproducible.
- RD is transformed with log10 for RMSE calculation.
- The two scores for your prediction of different well logs will be rank transformed, and averaged, e.g., if the **NMSE** of predictions for the values of RHOB, NPHI, and RD ranks 1st, and the **MAD** of assosiated depth shifts ranks 3rd, your final score would be  $$score = \frac{1+3}{2}=2 $$ The team with the minimum score win the competition
- Understanding and optimizing your predictions for this evaluation metric is paramount for this competition.


**Notes for submission:**
1. Only one user can register for the competition per team. 
2. __The user name has to be exactly the same as the team name__. If space is not allowed, please replace space with underscore '_'.
3. The submission file must be a zip file with three .csv files, which include your predictions of corrected/shfited well logs and their depth shift. **Note: the three csv files in your submission should have the same file names as the three files in the test set**, see name convention below.
- Zipped csv file name {team_name}\_submission\_{number of submission}.zip, which should contain (aligned_well_01.csv, aligned_well_02.csv, and aligned_well_03.csv)
- Each csv file should contain aligned well logs (NPHI_pred,RHOB_pred,RD_pred), which is calculated by shifting raw well logs (NPHI, RHOB, RD) according to reference depth (DEPT) and reference log (GR). The depth of raw well logs after depth shift (NPHI_dept_pred, RHOB_dept_pred,RD_dept_pred) also need to be included in the csv file for MAD calculation. Please refer to submission example "pdda_submission_1.zip" and don't change variable names.
4. The submission status might need a couple minutes to be updated, don't refresh the page too often.
5. The user needs to manually submit their best results to the leaderboard. Click "Participate", " Submit / View Results", click the "+" symbol in your submission. See the red circles in the attached figure.
6. Since Codalab doesn't support two scores or randk transformation, the score you see in the leaderboard is scaled for better display, the first 4 digits represent the scaled NMSE, and the last 4 digits represent the scaled MAD, the final ranking may be different from the leaderboard in Codalab after we collect all submissions. To see the actual scores for depth shift and well-log predictions, please Click "Participate", " Submit / View Results", click the "+" symbol in your submission, and Click "View scoring output log"
7. Please use version-control properly, as we need to validate your code and reproduce the results of the final submitted score in order to rank your team in the final scoreboard.  
8. Max submissions per day: 3
9. Max submissions total: 100



### <a>Timeline</a>

- __March 26, 2023__ - Registration starts. You must email Wen Pan (pdda_sig@spwla.org) with team information (team name, member names, affiliations, and emails) to get the submission URL.
- __March 31, 2023__ - Competition starts and data releases on github. 
- __August 26, 2023__ - Submission deadline. 
- __August 31, 2023__ - Announce winners.
- __September/October, 2023__ (tentative) - Award ceremony and presentations in the special session of the SPWLA Spring Topical Conference - Petrophysical Machine Learning.

All deadlines are at 11:59 PM UTC on the corresponding day unless otherwise noted. The competition organizers reserve the right to update the contest timeline if they deem it necessary.

### <a>Competition Rules</a>

1. Contestant can be an individual or a group with the maximum size of 4.
2. The contest focuses on data-driven methods, the use of additional data or petrophysical equations is not allowed.
3. Privately sharing code or data outside of teams is not permitted. However, it's okay to share code if made available to all participants on the competition Github repository via submitting issues or pull requests. 
4. A contestant will submit the predictions for both **corrected well log** and the **depth shift** of each type of well logs for each testing wells
5. A contestant will submit the source code and a brief report documenting the accuracy achieved in a few plots.
6. The judges will review the source code.
7. A leaderboard will be updating the rank of submissions from each team, but the score and rank is subject to scaling and change.
8. The contestant with the best quality source code and the best performance will be declared the winner for this competition.
    
### <a>Prize Policy</a>

- 1st Place - \$500  
- 2nd Place - \$400  
- 3rd Place - \$300   
- 4th Place - \$200   
- 5th Place - \$100   

Top 5 winning teams will be awarded with prizes(NOT in cash).

Note: The winners will additionally be required to provide a detailed description of their method in order to claim the prize (minimum of 2 pages double-column) by June 15, 2023.

Novel and practical algorithms will be recommended for a submission to the a SPWLA special issue by PDDA or a journal paper. 
<!-- #endregion -->

### <a>Contest Committee</a>
Wen Pan, Michael Ashby, Lei Fu, Yanxiang Yu, HyungJoo Lee, Chicheng Xu, Jaehyuk Lee 


