
### Check what the issue is:

# GAGTCGCGTA+AGCGAAGATT
# GAGTCGCTTC+AGCGAAGATT

#### Log in to bravolims:

ssh bravolims@abacus3


#### Edit the MM file

nano /lb/robot/research/processing/events/system/forcemm/HFMKGDSX3
# 1:1,1
# 2:2,1
# 3:1,1
# 4:1,1


#### Run the restart script to get the kill job commands, lane delete commands and event file copy commands:
#P.S Jobs preferably killed first below.

cd /home/bravolims/runprocessing
sh restartevent.sh "221012_A01433_0245_BHFMKGDSX3_MOHRun51"



#
# [bravolims@abacus3 runprocessing]$ sh restartevent.sh "221012_A01433_0245_BHFMKGDSX3_MOHRun51"
# 221012_A01433_0245_BHFMKGDSX3_MOHRun51
# 6 A01433_0245.listjobs#########] 503/503 (100%)
# ProcessLUID  ProjectLUID  ProjectName        ContainerLUID  ContainerName  Position  Index                                    LibraryLUID  LibraryProcess  ArtifactLUIDLibNorm  ArtifactNameLibNorm
# 24-401011    PAR1602      MOH-Q McGill_MUHC  27-76506       HFMKGDSX3      1:1       IDT_10nt_UDI_i7_185-IDT_10nt_UDI_i5_185  2-2415030    PCR-free        2-2512824            MoHQ-MU-8-19-FF1-1DT-IDT
# 24-401011    PAR1602      MOH-Q McGill_MUHC  27-76506       HFMKGDSX3      1:1       IDT_10nt_UDI_i7_184-IDT_10nt_UDI_i5_184  2-2415029    PCR-free        2-2512823            MoHQ-MU-8-18-FF1-1DT-IDT
# 24-401011    PAR1602      MOH-Q McGill_MUHC  27-76506       HFMKGDSX3      4:1       IDT_10nt_UDI_i7_191-IDT_10nt_UDI_i5_191  2-2415036    PCR-free        2-2512827            MoHQ-MU-8-22-FF1-1DT-IDT
# 24-401011    PAR1580      MOH-Q McGill_GCI   27-76506       HFMKGDSX3      4:1       IDT_10nt_UDI_i7_097-IDT_10nt_UDI_i5_097  2-2460566    PCR-free        2-2512828            MoHQ-GC-15-110-FT1-1DT-I
# 24-401011    PAR1602      MOH-Q McGill_MUHC  27-76506       HFMKGDSX3      3:1       IDT_10nt_UDI_i7_187-IDT_10nt_UDI_i5_187  2-2415032    PCR-free        2-2512825            MoHQ-MU-8-20-FF1-1DT-IDT
# 24-401011    PAR1580      MOH-Q McGill_GCI   27-76506       HFMKGDSX3      2:1       IDT_10nt_UDI_i7_099-IDT_10nt_UDI_i5_099  2-2460568    PCR-free        2-2512829            MoHQ-GC-15-109-FT1-1DT-I
# 24-401011    PAR1602      MOH-Q McGill_MUHC  27-76506       HFMKGDSX3      3:1       IDT_10nt_UDI_i7_189-IDT_10nt_UDI_i5_189  2-2415034    PCR-free        2-2512826            MoHQ-MU-8-21-FF1-1DT-IDT
# 24-401011    PAR1580      MOH-Q McGill_GCI   27-76506       HFMKGDSX3      2:1       IDT_10nt_UDI_i7_101-IDT_10nt_UDI_i5_101  2-2460570    PCR-free        2-2512830            MoHQ-GC-15-98-FT3-1DT-ID
# /lb/robot/research/processing/events/system/2022/2022-10-12-T11.26.29-valid/92-2530291_24-401011_samples.txt
# found event: 92-2530291_24-401011_samples.txt... New event -> start job.	Run end monitor process ID: 2392188 monitor job found: YES:
# (bravoli+ 2392188  6.7  0.0 115556  3300 pts/5    S+   Oct12 781:46 sh event_service.sh run)

  kill -9  2392188

# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq
# rm -r /lb/robot/research/processing/novaseq/2022/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq
#
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Unaligned.1
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Aligned.1
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/job_output/*/*.1.*done

rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Unaligned.2
rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Aligned.2
rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/job_output/*/*.2.*done


# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Unaligned.3
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Aligned.3
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/job_output/*/*.3.*done
#
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Unaligned.4
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/Aligned.4
# rm -r /nb/Research/processing/221012_A01433_0245_BHFMKGDSX3_MOHRun51-novaseq/job_output/*/*.4.*done

sh A01433_0245.listjobs
rm A01433_0245.listjobs


cp -v /lb/robot/research/processing/events/system/2022/2022-10-12-T11.26.29-valid/92-2530291_24-401011_samples.txt /lb/robot/research/processing/events/92-2530291_24-401011_samples_redo1.txt


