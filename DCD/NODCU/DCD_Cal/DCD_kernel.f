
       PROGRAM NODCU14
C     generate drainage, diversion, seepage, gw storage- for 142 islands, not nodes.
C     4 times of the original values in this code. 5/27/2016
C-----Suspected assigned leach water is much less than the actual value. Keep the original leach apl
C     and Drn files, 5 times of the original values in this code. 5/18/2016
C     If Oct-Nov have storms, the leach apl will reduce accordingly.  5/18/2016

C-----set groundwater rate as 0, remove the groundwater contribution   1/29/2016
c-----Add the groundwater daily output of 142 islands
c-----Change the monthly process into daily process and generate daily output  1/9/2013

c-----Added more output to file NODCU10-WY.21. Changes bracketed by CNM-10/04/95
c-----Added island diversion data without seepage (variable ID1)
c-----I'll use the data to compare DICU results to Report 4 AW data that
c-----did not include seepage.

c-----Modified program NODCU11                           NM (8/7/95)
c-----Modified to include Paul's new seepage estimates and irrigation efficiency by
c-----DOC subregion. See Chap 5 (Paul's work) in the 16th Annual Report.
c-----Seepage changes: A drained seepage component has been added.
c-----Irrigation efficiency: Varies by MWQI subregions.
c-----As a quick fix the drained seepage values will be added to both the div and returns.
c-----Both the seepage and irr eff. values are dependent on MWQI subregions.

c-----MWQI subregion   Drained seepage(AF/acre)     Irrigation efficiency
c-----High DOC region         0.095                    57%
c-----Midrange DOC            0.074                    67%
c-----Lower DOC region        0.013                    85%
c-----Changes are bracketed by C12

c-----Modified program NODCU10                          PNT (6/9/95)
c-----NODCU10 was used to process one water year at a time. The new version can
c-----process all the water years for which the hydrologic data is available
c-----in one run, resulting in considerable savings in CPU time. In addition
c-----a new ASCII file is created (junk.txt), which in turn can be used
c-----to convert this data into DSS format. To do this just type:

c-----         dssts <junk.txt

c-----This may take some time to process. After this conversion is completed,
c-----you may want to delete the file junk.txt, because it's fairly big.

c-----Modified program NODCU9                          NM (3/10/95)
c-----Modified to seperate seepage from irrigation diversions mainly
c-----for particle tracking modeling reasons. Total channel diversions will
c-----now be segregated and shown in two extra columns:irrigation diversions and seepage.
c-----Changes are bracketed by CNM
c-----ID1 is the irrigation diversion component of the total diversion from the channel
c-----ID2 is the seepage component of the total diversion from the channel
c-----Segregation of channel diversions was not calculated for ag den files and no
c-----changes were made to NODCU output files (NODCU10-WY.2?)


c-----Modified program NODCU8                          - NM (2/25/92)
c-----The main change is in the Ag Drain files. Parviz (PNT) needs
c-----a new format for the Ag Drain Model. Instead of P-ET, we want
c-----the change in soil moisture. This file soil.prc is created by a new
c-----version of dicu5 called dicu5.1. Changes can be noted by C@.
c-----Change hard coding of the number of nodes. Now, the number of nodes willbe read
c-----from the GEOM-NODES file.
c-----Changed DO LOOP maximum index for reading in DICU5.XX files

c-----PROGRAM NODCU8
c-----Modified NODCU7                                   - NM (5/7/91)

c----- Reads input files straight from dicu program. (eliminates pickwy.f)

c----- Asks whether drainage conc should be written in TDS or CL

c----- Output files have the water year as part of the filename

c----- The nodes written to the drn and div file are exactly the same as those
c        written to the drn quality file. ALL the nodes in the model are included
c        in the output file regardless of the drainage and diversion values. to do
c        this, node numbers (419 total) are read from a file called geom-nodes.
c
c----- Writes drn,div and precip-ET in Ag-Drain model input format.12 monthly files.
c        (All values in CFS). UNIT 30

C********************************************************************
C  MODIFIED NODCU6 :                                    - KG 11/06/89
C  *  RE-FORMAT OUTPUT FOR USE IN THE DELTA ISLAND MODEL
C  *  NEW PARAMETERS: WX,TAF2CFS,D2,DIV,DRN
C  *  NEW OUTPUT UNITS: 25,26,27
C  *  DETERMINE NET ATMOSPHERIC WATER EXCHANGE    ==> UNIT 25
C  *  DETERMINE NODAL DIVERSIONS BY ISLAND        ==> UNIT 26
C     INCLUDES SEEPAGE, APPLIED WATER, AND
C       APPLIED LEACH WATER
C  *  DETERMINE NODAL DRAINAGE RETURNS BY ISLAND  ==> UNIT 27
C     INCLUDES EXCESS APPLIED WATER AND
C       DRAINED LEACH WATER
C     DOES NOT INCLUDE PRECIP RUNOFF
C
C
C  MODIFIED NODCU5 :                                    - KG 09/13/89
C   * PREPARED FOR USING DRAINAGE TDS CONC FROM DWR BULLETIN 123
C     BY REGIONS (NORTH, SOUTHEAST, WEST)
C   * RE-ASSIGN ISLAND DRAINAGE TDS CONC :
C      ISL 102: NODE 274 SET FROM SOUTHEAST TO NORTH
C      ISL 103: NODES 253, 274, 276, 278 SET FROM SOUTHEAST TO NORTH
C      ISL 121: NODES 354, 355 SET FROM NORTH TO WEST
C      ISL 130: NODES 66 - 70, SET FROM WEST TO SOUTHEAST
C
C  MODIFIED NODCU3 :                                    - KG 06/21/89
C   * READ IN REPRESENTATIVE MONTHLY DRAINAGE TDS CONCENTRATION
C     FOR EACH 142 ISLAND AREAS (FILE IDRNTDS.DAT; UNIT 12; ARRAY DS)
C     POSSIBLE SOURCE: DWR BULLETIN 123
C   * ASSIGN NODAL DRAINAGE TDS CONCENTRATIONS
C     WEIGHTED BY ISLAND-TO-NODE DRAINAGE RETURN FLOW (ARRAY NDS)
C   * WRITE OUT NODAL DRAINAGE TDS CONCENTRATIONS
C     REPORT FORMAT ==> FILE NODCU5.24; UNIT 24
C     SEPARATE MONTHLY FILES : OCT-DTDS THRU SEP-DTDS (UNITS 51-62)
C
C  DETERMINATION OF IRRIGATION DIVERSIONS AND DRAINAGES FOR DELTA
C  ISLANDS AND ALLOCATION TO DWR/RMA DELTA MODEL NODES   - KG FEB. 88
C  INPUT FROM DICU5 RUN FOR CONSECUTIVE YEARS
C    -INCLUDING MONTHLY DELTA NET CHANNEL DEPLETIONS     - KG MAR. 88
C  RUNTIME SCREEN PRINT ADDED; INPUT FILE NAME CHANGE    - KG APR. 88
C
C  INPUT DATA INCLUDE:
C
C      1 - APPLIED WATER AND SEEPAGE (CONSUMPTIVE USE PROGRAM DICU4)
C      2 - LOWLANDS LEACHING REQUIREMENT ESTIMATES (G. SATO, 1981)
C      3 - NODAL ALLOCATION FACTORS (D. TAYLOR, 1987)
C      4 - IRRIGATION EFFICIENCY ESTIMATE (CALIBRATION PARAMETER)
C
C********************************************************************
C..
      CHARACTER*3 WMONTH
      CHARACTER*1 WYTYP
      CHARACTER*8 QUANTITY, tempconc
      CHARACTER*16 OUTFILE
      INTEGER*4 LPYR,YR,IFLAG,JIDFLAG,JDFLAG,CONC
      REAL*4 AW,S,IIE,LRA,LRD,FI,FD,ID,ID1,ID2,D,NID,NID1,NID2,ND
      REAL*4 LRAM, LRDM
      REAL*4 NIDSUB,NIDSUB1,NIDSUB2,NDSUB
      REAL*4 TAID,TAD,TMDID,TMDD,TADID,TADD,NTMDID,NTMDD,FISUM,FDSUM
      REAL*4 TCU,PREC,RO,CDNET,CDDNET,NCDNET
      REAL*4 NDS,TAF2CFS,WX,DIV,DIV1,DIV2,DRN,D2

C..
      INTEGER     MAXISL,MAXNODE
      PARAMETER   (MAXISL=168, MAXNODE=500)
      DIMENSION WMONTH(12),NODE(MAXNODE)
      DIMENSION DAYS_MONTH(12),DAYS_MONTH_LEAP(12)
C++++++++++++++++GROUNDWATER RATIOES++++++++++++++++++++++++++++++++
      DIMENSION GW_RATE(1922:2100)
C..
      COMMON/A/ YR,LPYR
      COMMON/B/ TADID,TADD,CDDNET
      COMMON/C/ TAID(MAXISL),TAD(MAXISL),TMDID(366),TMDD(366),
     .          CDNET(366),NTMDID(366),NTMDD(366),FISUM(MAXISL),
     &          FDSUM(MAXISL),IIE(MAXISL)
      COMMON/D/ S(MAXISL,366,1922:2100),AW(MAXISL,366,1922:2100),
     .          LRA(MAXISL,366),LRD(MAXISL,366),
     &          LRAM(MAXISL,12), LRDM(MAXISL, 12),
     &          FI(MAXNODE,MAXISL),
     &          FD(MAXNODE,MAXISL),TCU(MAXISL,366,1922:2100),
     .          PREC(MAXISL,366,1922:2100),RO(MAXISL,366),
     &          WNCU(MAXISL,366,1922:2100)
      COMMON/E/ ID(MAXISL,366),ID1(MAXISL,366),ID2(MAXISL,366),
     &          D(MAXISL,366),NID(MAXNODE,366),NID1(MAXNODE,366),
     .          NID2(MAXNODE,366),ND(MAXNODE,366),NIDSUB(MAXNODE,366),
     &          NIDSUB1(MAXNODE,366),NIDSUB2(MAXNODE,366),
     .          NDSUB(MAXNODE,366),NCDNET(366)
      COMMON/F/ DS(MAXISL,366),NDS(MAXNODE,366)
      COMMON/G/ WX(MAXISL,366),D2(MAXISL,366),DIV(MAXISL,366),
     &          DIV1(MAXISL,366),DIV2(MAXISL,366),DRN(MAXISL,366)
      COMMON/H/ GWAMOUNT1(MAXISL,366),GWAMOUNT2(MAXISL,366),
     &          GWINP(MAXISL,6),GWF(MAXISL,366,1922:2100),
     &          DRN_Y(MAXISL,366,1922:2100),DIV_Y(MAXISL,366,1922:2100),
     &          SPG_Y(MAXISL,366,1922:2100),RO_Y(MAXISL,366,1922:2100),
     &        DSPG_Y(MAXISL,366,1922:2100),LRD_Y(MAXISL,366,1922:2100),
     &        LRA_Y(MAXISL,366,1922:2100),tempD_Y(MAXISL,366,1922:2100),
     &          GWA2_Y(MAXISL,366,1922:2100)
C..
      DATA WMONTH/'OCT','NOV','DEC','JAN','FEB','MAR','APR','MAY',
     .            'JUN','JUL','AUG','SEP'/
      DATA DAYS_MONTH/31,61,92,123,151,182,212,243,273,304,335,365/
      DATA DAYS_MONTH_LEAP/31,61,92,123,152,183,213,244,274,305,336,366/

      CHARACTER*1 JUNK, ANSWER
	CHARACTER*4 tempyear
	CHARACTER*80 DICU5_14,DICU5_17,DICU5_12,DICU5_27,DICU5_30,DIVFCTR_RMA,
     &	           DRNFCTR_RMA,LEACHAPL_DAT, GROUNDWATER
	CHARACTER*80 LEACHDRN_DAT,IDRNTDS_DAT,DRNCL_123,GEOM_NODES,IRREFF_DAT,
     &             subarea_info

      REAL QDRN(MAXNODE,1922:2100,366),QIRR(MAXNODE,1922:2100,366),
     &     QSEEP(MAXNODE,1922:2100,366)
      REAL WYR(1922:2100),TEMPYR
C***************************MONTHS TO DAYS************************************
      LOGICAL NDFLG(MAXNODE)

C12A
      CHARACTER*1 docregion(maxisl)
      CHARACTER*100 junk1
	CHARACTER*120 FLE
      INTEGER uplow(maxisl),area(maxisl),drnseep(maxisl)
C12B

C..
      WRITE(*,'(48H*** PROGRAM NODCU12---VERSION January 2013*** //)')

      TAF2CFS = 0.50417
C     actually TAF2CFS here is AF2CFS (43560/24/3600)
      DO N=1,MAXNODE
         NDFLG(N)=.FALSE.
      ENDDO

      IYR1=1922
      IYR2=2100
C++++++++++++++++GROUNDWATER RATIOES++++++++++++++++++++++++++++++++
      DO YR=IYR1,IYR2
         GW_RATE(YR)=0.0
      ENDDO

 901  FORMAT(' THE CURRENT VERSION OF NODCU PROGRAM ASSUMES THAT THE'/
     &     ' DATABASE CONTAINS DATA FROM  YEAR: ',I4,'  TO  ',I4//
     &     ' IS THIS ASSUMPTION CORRECT? (Y/N): ',$)
ccccccccccccccccccccccccccccccccccc      READ(*,801)ANSWER
	CALL GETENV('years_ok',ANSWER)

      IF(ANSWER.EQ.'N' .OR. ANSWER.EQ.'n')THEN
         WRITE(*,902)
 902     FORMAT(' ENTER BEGINNING AND ENDING YEAR yyyy yyyy: ',$)
cccccccccccccccccccccccccccccccccc         READ*,IYR1,IYR2
 	   CALL GETENV('begwy',tempyear)
       READ(tempyear,'(I4)') IYR1
       CALL GETENV('endwy',tempyear)
       READ(tempyear,'(I4)') IYR2
	   CALL GETENV('leachscale',tempconc)
	   READ(tempconc,'(I8)') leachscale
      ENDIF

      DO I=1,MAXISL
         TAID(I)=0.0
         TAD(I)=0.0
         FISUM(I)=0.0
         FDSUM(I)=0.0
         IIE(I)=0.0
         DO N=1,MAXNODE
            FI(N,I)=0.0
            FD(N,I)=0.0
         ENDDO
      ENDDO

      DO M=1,366
         TMDID(M)=0.0
         TMDD(M)=0.0
         CDNET(M)=0.0
         NTMDID(M)=0.0
         NTMDD(M)=0.0
         NCDNET(M)=0.0
         DO I=1,MAXISL
            DIV(I,M)=0.0
            DIV1(I,M)=0.0
            DIV2(I,M)=0.0
            DRN(I,M)=0.0
            LRA(I,M)=0.0
            LRD(I,M)=0.0
            RO(I,M)=0.0
            WX(I,M)=0.0
            ID(I,M)=0.0
            ID1(I,M)=0.0
            ID2(I,M)=0.0
            D(I,M)=0.0
            D2(I,M)=0.0
            DS(I,M)=0.0
            DO YR=IYR1,IYR2
               S(I,M,YR)=0.0
               AW(I,M,YR)=0.0
               TCU(I,M,YR)=0.0
               PREC(I,M,YR)=0.0
               WNCU(I,M,YR)=0.0
            ENDDO
         ENDDO

         DO N=1,MAXNODE
            NIDSUB(N,M)=0.0
            NIDSUB1(N,M)=0.0
            NIDSUB2(N,M)=0.0
            NDSUB(N,M)=0.0
            NID(N,M)=0.0
            NID1(N,M)=0.0
            NID2(N,M)=0.0
            ND(N,M)=0.0
            NDS(N,M)=0.0
         ENDDO
      ENDDO
Ccccccccccccccccccccccccccccccccccccccccccccccccccc      READ(*,*) CONC
	CALL GETENV('datatype',tempconc)
      READ(tempconc,'(I4)') CONC

      WRITE(*,951)
 951  FORMAT(/' DO YOU WANT TO CREATE AN ASCII FILE FOR LOADING DATA
     &   IN DSS?(Y/N): ',$)
CccccccccccccccccccccccccccccccCcccccccccccccccccccccc      READ(*,801)ANSWER
	CALL GETENV('ascii',ANSWER)
      write(*,*) "  "
      IF(ANSWER.EQ.'Y' .OR. ANSWER.EQ.'y')THEN
Ccccccccccccccccccccccccccccccccccccccccccccccccccc         READ(*,929)FLE
	   CALL GETENV('dssfile',FLE)
	   write(*,*) FLE
c 929     FORMAT(A120)
      ENDIF



C-----1. READS INPUT FILES FROM THE NEW PICKWY PROGRAM (YEAR IS PART OF FILENAME)
	CALL GETENV('DICU5_14',DICU5_14)
	CALL GETENV('DICU5_17',DICU5_17)
	CALL GETENV('IRREFF_DAT',IRREFF_DAT)
	CALL GETENV('LEACHAPL_DAT',LEACHAPL_DAT)
	CALL GETENV('LEACHDRN_DAT',LEACHDRN_DAT)
	CALL GETENV('DICU5_12',DICU5_12)
	CALL GETENV('DICU5_27',DICU5_27)
	CALL GETENV('DICU5_30',DICU5_30)
	CALL GETENV('GW_RATES_TXT',GROUNDWATER)

	OPEN(1,FILE=DICU5_14,FORM='FORMATTED',STATUS='OLD')
      OPEN(2,FILE=DICU5_17,FORM='FORMATTED',STATUS='OLD')
      OPEN(3,FILE=IRREFF_DAT,FORM='FORMATTED',STATUS='OLD')
      OPEN(8,FILE=LEACHAPL_DAT,FORM='FORMATTED',STATUS='OLD')
      OPEN(9,FILE=LEACHDRN_DAT,FORM='FORMATTED',STATUS='OLD')
      OPEN(10,FILE=DICU5_12,FORM='FORMATTED',STATUS='OLD')
      OPEN(11,FILE=DICU5_27,FORM='FORMATTED',STATUS='OLD')
      OPEN(85,FILE=DICU5_30,FORM='FORMATTED',STATUS='OLD')
C++++++++++++++++GROUNDWATER RATIOES++++++++++++++++++++++++++++++++
      OPEN(15,FILE=GROUNDWATER,STATUS='OLD')
      open(17,file='WYTYPES',form='formatted',status='old')
      open(56,file='gwbyisl.txt',status='unknown')
      open(58,file='drn_wo_ro_isl.txt',status='unknown')
      open(60,file='div_wo_spg_isl.txt',status='unknown')
      open(62,file='spgisl.txt',status='unknown')
      open(64,file='roisl.txt',status='unknown')
      open(66,file='LRDisl.txt',status='unknown')
      open(68,file='LRAisl.txt',status='unknown')
      open(70,file='tempdrnsisl.txt',status='unknown')
      open(72,file='GWAmount2isl.txt',status='unknown')

      DO YR = IYR1,IYR2
        READ(17,*) WYTYP,TEMPYR
        IF(WYTYP.EQ.'C'.OR.WYTYP.EQ.'D') THEN
            WYR(YR) = 1
        ELSE
            WYR(YR) = 0
        ENDIF
      ENDDO

C*    MODIFY PROGRAM TO OUTPUT EITHER CL OR TDS CONCENTRATIONS --NM(5/7/91)
      IF (CONC.EQ.1) THEN
	CALL GETENV('IDRNTDS_DAT',IDRNTDS_DAT)
	write(*,*) "what is the IDRNTDS =", IDRNTDS_DAT
      OPEN(12,FILE=IDRNTDS_DAT,FORM='FORMATTED',STATUS='OLD')
      ENDIF

      IF (CONC.EQ.2) THEN
	CALL GETENV('DRNCL_123',DRNCL_123)
      OPEN(12,FILE=DRNCL_123,FORM='FORMATTED',STATUS='OLD')
      ENDIF


C-----write the island diverion, drainage and seepages
      OPEN(80, FILE="ISL_DIV.txt")
      OPEN(81, FILE="ISL_DRN.txt")
      OPEN(82, FILE="ISL_SPG.txt")

C++++++++++++++++GROUNDWATER RATIOES++++++++++++++++++++++++++++++++
      DO YR = IYR1,IYR2
         READ(15,*) IYR, GW_RATE(YR)
      ENDDO

      DO ISL=1,MAXISL
         READ(2,910)JUNK
         READ(1,910)JUNK
         READ(10,910)JUNK
         READ(11,910)JUNK
         READ(85,910)JUNK
 910     FORMAT(/A1)

         DO YR=IYR1,IYR2
C********************************************************
            IF (MOD(YR,4).EQ.0) THEN
                DAYS_YR = 366
            ELSE
                DAYS_YR = 365
            ENDIF
C*********LEAP YEAR OR NOT, DAYS OF THE YEAR************

            READ(2,980) ISLDICU,IYR,(AW(ISLDICU,M,YR),M=1,DAYS_YR)
            IF (IYR.NE.YR.OR. ISL.NE.ISLDICU) THEN
               WRITE(*,990)
 990           FORMAT(' ERROR... Changes in the format of DICU output
     &          files'/ '          are detected. NODCU program must
     &          be modified'/  '          to accomodate the changes')
               STOP
            ENDIF
            READ(1,980) ISLDICU,IYR,(S(ISLDICU,M,YR),M=1,DAYS_YR)
            READ(10,980)ISLDICU,IYR,(PREC(ISLDICU,M,YR),M=1,DAYS_YR)
            READ(11,980) ISLDICU,IYR,(TCU(ISLDICU,M,YR),M=1,DAYS_YR)
            READ(85,980)ISLDICU,IYR,(WNCU(ISLDICU,M,YR),M=1,DAYS_YR)

 980		  FORMAT(I3,4X,I4,366F8.0)
         ENDDO
         IF(MOD(ISL,10).EQ.0 .OR. ISL.EQ.MAXISL)THEN
            WRITE(*,991)ISL
 991        FORMAT(' FINISHED READING DATA FOR SUBAREA: ',I4)
         ENDIF
      ENDDO

      READ(3,'(f4.2)')(IIE(I),I=1,MAXISL)
      READ(8,15)((LRAM(I,M),M=1,12),I=1,MAXISL)
      READ(9,15)((LRDM(I,M),M=1,12),I=1,MAXISL)
      do i = 1, 12
        write(*,*)"LRA ", i," = ", LRAM(1,i)
      enddo
C12A---Read in data to calculate new seepage component (drained seepage)

	CALL GETENV('subarea_info',subarea_info)
      open(100,file=subarea_info,form='formatted',status='old')

      do j=1,6
         read(100,'(a100)') junk1
      enddo

      do i=1,maxisl
         read(100,'(57X,I1,5X,I10,8X,A1)') uplow(i),area(i),docregion(i)
      end do

      close(100)

C---Calculte drained seepage for each subarea

      do 100 i=1,maxisl

         if (uplow(i).eq.2.or.i.eq.133.or.i.eq.134.or.i.eq.135.or.i.eq.
     &     136.or.i.eq.137.or.i.eq.140.or.i.eq.141.or.i.eq.142) then
            drnseep(i)=0
            go to 100
         endif
         if (uplow(i).eq.1) then
            if (docregion(i).eq."L") drnseep(i)=0.013*area(i)
            if (docregion(i).eq."M") drnseep(i)=0.074*area(i)
            if (docregion(i).eq."H") drnseep(i)=0.095*area(i)
            go to 100
         endif
         write(*,*) " Aghhhh! Subarea is neither in the Delta Lowlands
     &      or Uplands"
         stop
 100   continue

C--------Add drn seepage component to seepage calculated by DICU program amd reported in DICU5.14

         do yr=iyr1,iyr2
C********************************************************
            IF(MOD(YR,4).EQ.0) LPYR=1
            IF (LPYR.EQ.1) THEN
                DAYS_YR = 366
                DO I = 1,MAXISL
                    DO M= 1,DAYS_YR
                        DO MM = 1,12
                            IF(MM.EQ.1) THEN
                                IF(M.LE.DAYS_MONTH_LEAP(MM)) THEN
                                s(i,m,yr)=s(i,m,yr)+drnseep(i)
     &                           /DAYS_MONTH_LEAP(MM)
                                ENDIF
                            ELSE
                                IF(M.GT.DAYS_MONTH_LEAP(MM-1).AND.
     &                          M.LE.DAYS_MONTH_LEAP(MM)) THEN
                                s(i,m,yr)=s(i,m,yr)+drnseep(i)
     &                           /(DAYS_MONTH_LEAP(MM)
     &                           -DAYS_MONTH_LEAP(MM-1))
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
             ELSE
                DAYS_YR = 365
                DO I = 1,MAXISL
                    DO M= 1,DAYS_YR
                        DO MM = 1,12
                            IF (MM.EQ.1) THEN
                                IF(M.LE.DAYS_MONTH(MM)) THEN
                                    s(i,m,yr)=s(i,m,yr)+drnseep(i)
     &                               /DAYS_MONTH(MM)
                               ENDIF
                            ELSE
                                IF(M.GT.DAYS_MONTH(MM-1).AND.
     &                             M.LE.DAYS_MONTH(MM)) THEN
                                    s(i,m,yr)=s(i,m,yr)+drnseep(i)/
     &                               (DAYS_MONTH(MM)-DAYS_MONTH(MM-1))
                                ENDIF
                           ENDIF
                        ENDDO
                    ENDDO
                ENDDO
             ENDIF
         enddo


C-----READ THE LIST OF DWRDSM MODEL NODES

	CALL GETENV('GEOM_NODES',GEOM_NODES)
      OPEN(50,FILE=GEOM_NODES,FORM='FORMATTED',STATUS="OLD")
      READ(50,*) NUMBER
      DO 805 I=1,NUMBER
         READ (50,*) NODE(I)
 805  CONTINUE

    8 CONTINUE
      READ(12,19)
      READ(12,13) ((DS(I,M),M=1,12),I=1,MAXISL)

C*
      WRITE(*,*) " FINISHED READING AND WRITING INPUT DATA"

C..
C  NODAL ALLOCATION FACTOR CHECK
C..
      DO 150 I=1,MAXISL
      DO 125 N=1,MAXNODE
      FISUM(I)=FISUM(I)+FI(N,I)
      FDSUM(I)=FDSUM(I)+FD(N,I)
  125 CONTINUE
  150 CONTINUE
C..
C..
C  DETERMINE PRECIPITATION RUNOFF IN ACRE-FT
C
C  DETERMINE NET ATMOSPHERIC WATER EXCHANGE IN ACRE-FT
C     POSITIVE ==> NET WATER ADDED TO ISLAND FROM ATMOSPHERE
C     NEGATIVE ==> NET WATER REMOVED FROM ISLAND TO ATMOSPHERE
C..
      DO YR=IYR1,IYR2
         TADID=0.
         TADD=0.
         CDDNET=0.
         print*,' processing yr:',YR
         LPYR=0
         RealYR=float(YR)
         IF(MOD(RealYR,4.).EQ.0) LPYR=1
         IF (LPYR.EQ.1) THEN
            DAYS_YR = 366
            DO 205 I = 1,MAXISL
                DO 204 M= 1,DAYS_YR
                    DO 203 MM = 1,12
                        IF(MM.EQ.1) THEN
                            IF(M.LE.DAYS_MONTH_LEAP(MM)) THEN
                            LRA(I,M) = LRAM(I,MM)/DAYS_MONTH_LEAP(MM)
                            LRD(I,M) = LRDM(I,MM)/DAYS_MONTH_LEAP(MM)
                            ENDIF
                        ELSE
                            IF(M.GT.DAYS_MONTH_LEAP(MM-1).AND.
     &                      M.LE.DAYS_MONTH_LEAP(MM)) THEN
                            LRA(I,M) = LRAM(I,MM)/(DAYS_MONTH_LEAP(MM)
     &                       -DAYS_MONTH_LEAP(MM-1))
                            LRD(I,M) = LRDM(I,MM)/(DAYS_MONTH_LEAP(MM)
     &                       -DAYS_MONTH_LEAP(MM-1))
                            ENDIF
                        ENDIF
 203                CONTINUE
 204        CONTINUE
 205        CONTINUE
         ELSE
            DAYS_YR = 365
            DO 208 I = 1,MAXISL
                DO 206 M= 1,DAYS_YR
                    DO 207 MM = 1,12
                        IF (MM.EQ.1) THEN
                            IF(M.LE.DAYS_MONTH(MM)) THEN
                                LRA(I,M) = LRAM(I,MM)/DAYS_MONTH(MM)
                                LRD(I,M) = LRDM(I,MM)/DAYS_MONTH(MM)
                            ENDIF
                        ELSE
                            IF(M.GT.DAYS_MONTH(MM-1).AND.
     &                       M.LE.DAYS_MONTH(MM)) THEN
                                LRA(I,M) = LRAM(I,MM)/(DAYS_MONTH(MM)
     &                           -DAYS_MONTH(MM-1))
                                LRD(I,M) = LRDM(I,MM)/(DAYS_MONTH(MM)
     &                          -DAYS_MONTH(MM-1))
                            ENDIF
                       ENDIF
 207            CONTINUE
 206        CONTINUE
 208        CONTINUE
         ENDIF

         DO 300 I=1,MAXISL
            DO 200 M=1,DAYS_YR
                WX(I,M)=PREC(I,M,YR)
                RO(I,M)=PREC(I,M,YR)
               IF(RO(I,M).LT.0) THEN
                  RO(I,M)=0.
               ENDIF
 200        CONTINUE
 300     CONTINUE
C..
         DO 99 I=20,27
            IF(YR.NE.IYR1) CLOSE(I)
			PRINT *,YR,IYR1,I
 99      CONTINUE


         DO 500 I=1,MAXISL
            DO 531 IG = 1, DAYS_YR
             if (uplow(i).eq.1) then
                if (docregion(i).eq."L") gwrate=0.35 !GW_RATE(YR)
                if (docregion(i).eq."M") gwrate=0.30 !0.4
                if (docregion(i).eq."H") gwrate=0.25 !0.4
             else
                gwrate = GW_RATE(YR)
             endif
               GWAMOUNT1(I,IG) = gwrate*AW(I,IG,YR)/IIE(I)
               GWAMOUNT2(I,IG) = gwrate*S(I,IG,YR)

               GWF(I,IG,YR)= GWAMOUNT1(I,IG)+GWAMOUNT2(I,IG)

 531        CONTINUE


            leachreduced = 0

            DO 400 M=1,DAYS_YR

               RO(I,M) = 0.75*RO(I,M)

               if (YR.eq.1977.and.M.ge.274.and.M.le.303) then
                  if (I.eq.132) then
                    AW(I,M,YR) = AW(I,M,YR)+ 430/30.
                  elseif (I.eq.61) then
                    AW(I,M,YR) = AW(I,M,YR)+ 700/30.
                  elseif (I.eq.57) then
                    AW(I,M,YR) = AW(I,M,YR) + 50/30.
                  elseif (I.eq.130) then
                    AW(I,M,YR) = AW(I,M,YR)+ 98/30.
                  endif
               endif

C   drained seepage and groundwater to lower the groundwater table
                tempe = (1.-IIE(I))/IIE(I)
               IF(MOD(YR,4).EQ.0) LPYR=1
               IF (LPYR.EQ.1) THEN

                      IF(M.GT.152) THEN
                          DO MM = 6,12
                            IF(M.GT.DAYS_MONTH_LEAP(MM-1).AND.
     &                        M.LE.DAYS_MONTH_LEAP(MM)) THEN
                              tempdrns = drnseep(i)/(DAYS_MONTH_LEAP(MM)
     &                          -DAYS_MONTH_LEAP(MM-1))
                            ENDIF
C                            IF(tempdrns.LT.((AW(I,M,YR)*tempe)+LRD(I,M)
C     &                       +RO(I,M))) THEN
C                                tempdrns = 0
C                            ENDIF
                          ENDDO
                      ELSE
                          IF(M.LE.DAYS_MONTH_LEAP(1)) THEN
                               tempdrns = drnseep(i)/DAYS_MONTH_LEAP(1)
C                                tempdrns = 0
                          ELSE
                            DO MM = 2,5
                              IF(M.GT.DAYS_MONTH_LEAP(MM-1).AND.
     &                        M.LE.DAYS_MONTH_LEAP(MM)) THEN
                              tempdrns = drnseep(i)/(DAYS_MONTH_LEAP(MM)
     &                          -DAYS_MONTH_LEAP(MM-1))
C                               tempdrns = 0
                              ENDIF
                            ENDDO
                          ENDIF
                      ENDIF


               ELSE
                  IF (M.GT.151) THEN
                       DO MM = 6,12
                            IF(M.GT.DAYS_MONTH(MM-1).AND.
     &                        M.LE.DAYS_MONTH(MM)) THEN
                              tempdrns = drnseep(i)/(DAYS_MONTH(MM)
     &                          -DAYS_MONTH(MM-1))
                            ENDIF
C                            IF(tempdrns.LT.((AW(I,M,YR)*tempe)+LRD(I,M)
C     &                       +RO(I,M))) THEN
C                                tempdrns = 0
C                            ENDIF
                      ENDDO
                  ELSE
                      IF(M.LE.DAYS_MONTH(1)) THEN
                          tempdrns = drnseep(i)/DAYS_MONTH(1)
C                          tempdrns = 0
                      ELSE
                          DO MM = 2,5
                            IF(M.GT.DAYS_MONTH(MM-1).AND.
     &                          M.LE.DAYS_MONTH(MM)) THEN
C                                tempdrns = 0
                                tempdrns = drnseep(i)/
     &                              (DAYS_MONTH(MM)-DAYS_MONTH(MM-1))
                            ENDIF
                          ENDDO
                      ENDIF
                  ENDIF
                 ENDIF



C +++++++++++++Leach water apl and drn reassumed double the originals and delay the apl time

               LRA(I,M) = leachscale*LRA(I,M)
               LRD(I,M) = 1.0*LRD(I,M)
               if (LRA(I,M).gt.RO(I,M)) then
                    if(leachreduced.gt.0.0001) then
                        if(leachreduced.gt.LRA(I,M)) then
                            LRA(I,M) = 0.0
                            leachreduced = leachreduced - LRA(I,M)
                        else
                            LRA(I,M) = LRA(I,M) - leachreduced
                            leachreduced = 0.0
                        endif
                    else
                        LRA(I,M) = LRA(I,M)- RO(I,M)
                        leachreduced = leachreduced+RO(I,M)
                    endif
               else
                    LRA(I,M) = 0.
                    leachreduced = leachreduced+LRA(I,M)
               endif
               if(LRD(I,M).gt.0.and.leachreduced.gt.0.0001)then
                    if (M.gt.92.and.M.le.123) then
                        LRD(I,M) = LRD(I,M)- leachreduced*0.56/31
                    else if (M.gt.123.and.M.le.151) then
                        LRD(I,M) = LRD(I,M)- leachreduced*0.29/28
                    else if(M.gt.151.and.M.le.182) then
                        LRD(I,M) = LRD(I,M)- leachreduced*0.14/31
                    else
                        LRD(I,M) = LRD(I,M)- leachreduced*0.01/31
                    endif
                    if (LRD(I,M).lt.0.0) LRD(I,M)= 0.0
               endif

               ID(I,M)=S(I,M,YR)-GWAMOUNT2(I,M)+LRA(I,M)+
     &                +AW(I,M,YR)/IIE(I)-GWAMOUNT1(I,M)+WNCU(I,M,YR)
CNM1-----------
              ID1(I,M) = AW(I,M,YR)/IIE(I)-GWAMOUNT1(I,M)+WNCU(I,M,YR)
     &                +LRA(I,M)
C               ID1(I,M)= ((1-gwrate)*AW(I,M,YR)/IIE(I)) + LRA(I,M)
C     &          + WNCU(I,M,YR)
              ID2(I,M) = S(I,M,YR)-GWAMOUNT2(I,M)
              IF(ID(I,M).lt.0) ID(I,M) = 0
              IF(ID1(I,M).lt.0) ID1(I,M)= 0.0
              IF(ID2(I,M).lt.0) ID2(I,M)=0.0
CNM2-----------

C12A------------Add the drained seepage component to the drainage also

               tempe = (1.-IIE(I))/IIE(I)

               D(I,M)=( AW(I,M,YR) * tempe ) + LRD(I,M)
     &                 + RO(I,M) + tempdrns
               D2(I,M)=( AW(I,M,YR) * tempe ) + LRD(I,M)
     &                 + tempdrns

C12B
                DRN_Y(I,M,YR) = D2(I,M)
                LRD_Y(I,M,YR) = LRD(I,M)
                DIV_Y(I,M,YR) = ID1(I,M)
                LRA_Y(I,M,YR) = LRA(I,M)
                SPG_Y(I,M,YR)= S(I,M,YR)-GWAMOUNT2(I,M)
                tempD_Y(I,M,YR) = tempdrns
                RO_Y(I,M,YR) = RO(I,M)
                GWA2_Y(I,M,YR) = GWAMOUNT2(I,M)
 400        CONTINUE
 500     CONTINUE
C..
C  SUMMATIONS BY ISLAND AND DELTA OVER ALL MONTHS
C..
         DO 550 I=1,MAXISL
            TAID(I)=0.
            TAD(I)=0.
            DO 525 M=1,DAYS_YR
               TAID(I)=TAID(I) + ID(I,M)
               TAD(I)=TAD(I) + D(I,M)
 525        CONTINUE
            TADID=TADID + TAID(I)
            TADD=TADD + TAD(I)
 550     CONTINUE
         CDDNET=TADID-TADD
C..
C  SUMMATIONS BY MONTH OVER ALL ISLANDS
C..
         DO 600 M=1,DAYS_YR
            TMDID(M)=0.
            TMDD(M)=0.
            DO 575 I=1,MAXISL
               TMDID(M)=TMDID(M) + ID(I,M)
               TMDD(M)=TMDD(M) + D(I,M)
 575        CONTINUE
            CDNET(M)=TMDID(M)-TMDD(M)
 600     CONTINUE
C..
C  WRITE OUT ISLAND RESULTS IN ACRE-FEET
C..
C  add the island outputs
         WRITE(80,*)YR
         WRITE(80,11)(I,(ID1(I,M),M=1,DAYS_YR),I=1,MAXISL)
         WRITE(81,*)YR
         WRITE(81,11)(I,(D(I,M),M=1,DAYS_YR),I=1,MAXISL)
         WRITE(82,*)YR
         WRITE(82,11)(I,(ID2(I,M),M=1,DAYS_YR),I=1,MAXISL)
         write(*,*) YR
CNM-10/04/95

C..
C  ALLOCATE TO DWR/DSM DELTA MODEL NODES
C  USING ISLAND-TO-NODE ALLOCATION FACTORS
C
C  DETERMINE WEIGHTED NODAL DRAINAGE TDS CONCENTRATIONS - KG 06/20/89
C  USING ISLAND DRAINAGE TDS CONCENTRATIONS IN ARRAY DS
C  AND ISLAND VOLUME CONTRIBUTIONS TO EACH NODE (NDSUB);
C  LOAD ARRAY NDS
C..

         DO 750 N=1,MAXNODE
            DO 700 M=1,DAYS_YR
               NID(N,M)=0.
CNM1-----------
               NID1(N,M)=0.
               NID2(N,M)=0.
CNM2-----------
               ND(N,M)=0.
               NDS(N,M)=0.
               DO 650 I=1,MAXISL
                  NIDSUB(N,M)=FI(N,I)/100. * ID(I,M)

CNM1--------------
                  NIDSUB1(N,M)=FI(N,I)/100. * ID1(I,M)
                  NIDSUB2(N,M)=FI(N,I)/100. * ID2(I,M)
CNM2--------------


                  NDSUB(N,M)=FD(N,I)/100. * D(I,M)
C..---------------
                  NID(N,M)=NID(N,M) + NIDSUB(N,M)
CNM1--------------

                  NID1(N,M)=NID1(N,M) + NIDSUB1(N,M)
                  NID2(N,M)=NID2(N,M) + NIDSUB2(N,M)
CNM2--------------

                  ND(N,M)=ND(N,M) + NDSUB(N,M)
C...
C  REASSIGN SELECTED ISLAND DRAINAGE CONC           - KG 09/14/89
C  CONSISTANT WITH BULLETIN 123
C...
                  IF(I.EQ.102 .AND. N.EQ.274) THEN
                     DSZ=DS(9,M)
                     NDS(N,M)=NDS(N,M) + (NDSUB(N,M)*DSZ)
                  ELSEIF((I.EQ.103) .AND. (N.EQ.253 .OR. N.EQ.274 .OR.
     >                    N.EQ.276 .OR.N.EQ.278)) THEN
                     DSZ=DS(9,M)
                     NDS(N,M)=NDS(N,M) + (NDSUB(N,M)*DSZ)
                  ELSEIF((I.EQ.121) .AND. (N.EQ.354 .OR. N.EQ.355)) THEN
                     DSZ=DS(8,M)
                     NDS(N,M)=NDS(N,M) + (NDSUB(N,M)*DSZ)
                  ELSEIF((I.EQ.130) .AND. (N.GE.66 .AND. N.LE.70)) THEN
                     DSZ=DS(31,M)
                     NDS(N,M)=NDS(N,M) + (NDSUB(N,M)*DSZ)
                  ELSE
C...
C...
                     NDS(N,M)=NDS(N,M) + (NDSUB(N,M)*DS(I,M))
                  ENDIF
 650           CONTINUE
               IF(ND(N,M).GT.0.) NDS(N,M)=NDS(N,M)/ND(N,M)
 700        CONTINUE
 750     CONTINUE
C..
C  SUM UP OVER ALL NODES FOR EACH MONTH
C..
         DO 785 M=1,DAYS_YR
            NTMDID(M)=0.
            NTMDD(M)=0.
            DO 775 N=1,MAXNODE
               NTMDID(M)=NTMDID(M) + NID(N,M)
               NTMDD(M)=NTMDD(M) + ND(N,M)
 775        CONTINUE
            NCDNET(M)=NTMDID(M)-NTMDD(M)
 785     CONTINUE

C  CONVERT FROM ACRE-FEET TO MEAN MONTHLY CFS
C  ALLOWING FOR LEAP YEARS
C..
         DO 900 N=1,MAXNODE
            DO 800 M=1,DAYS_YR
               NID(N,M)=NID(N,M) * TAF2CFS
               NID1(N,M)=NID1(N,M) * TAF2CFS
               NID2(N,M)=NID2(N,M) * TAF2CFS
               ND(N,M)=ND(N,M) * TAF2CFS
 800        CONTINUE
 900     CONTINUE
C..
         DO 950 M=1,DAYS_YR
            NTMDID(M)=NTMDID(M) * TAF2CFS
            NTMDD(M)=NTMDD(M) * TAF2CFS
            NCDNET(M)=NCDNET(M) * TAF2CFS
 950     CONTINUE

         DO 955 N=1,MAXNODE
            IFLAG=0
            DO 953 M=1,DAYS_YR
               IF(NDS(N,M).NE.0.) IFLAG=IFLAG+1
 953        CONTINUE
CC            IF(IFLAG.NE.0) WRITE(24,1049) N,(NDS(N,M),M=1,DAYS_YR)
 955     CONTINUE
C...
C...
C  DETERMINE NODAL DIVERSIONS AND DRAINAGE RETURNS BY ISLAND AND NODE
C  FOR THE DELTA ISLAND MODEL
C
C  WRITE NON ZERO DIVERSIONS BY ISLAND       TO UNIT 26
C  WRITE NON ZERO DRAINAGE RETURNS BY ISLAND TO UNIT 27
C...
C         TAF2CFS(5)=0.018
C         IF(LPYR.EQ.1) TAF2CFS(5)=0.01738

C...
         DO 780 I=1,MAXISL
            DO 776 N=1,MAXNODE
               JIDFLAG=0
               JDFLAG=0
               DO 765 M=1,366
C-----------------DIV(M)=0.0
C-----------------DRN(M)=0.0
 765           CONTINUE
               DO 770 M=1,DAYS_YR
                  DIV(I,M)=FI(N,I)/100. * ID(I,M) * TAF2CFS

CNM1--------------
                  DIV1(I,M)=FI(N,I)/100. * ID1(I,M) * TAF2CFS
                  DIV2(I,M)=FI(N,I)/100. * ID2(I,M) * TAF2CFS
CNM2--------------


                  DRN(I,M)=FD(N,I)/100. * D2(I,M) * TAF2CFS
                  IF(DIV(I,M).NE.0.) JIDFLAG=1
CNM1--------------
                  IF(DIV1(I,M).NE.0.) JIDFLAG=1
                  IF(DIV2(I,M).NE.0.) JIDFLAG=1
CNM2--------------
                  IF(DRN(I,M).NE.0.) JDFLAG=1
 770           CONTINUE
 776        CONTINUE
 780     CONTINUE

C-----4. FILES WRITTEN FOR AG-DRAIN (DELTA ISLAND) MODEL.

C..
C  RESULTS IN DWR/DSM DELTA MODEL INPUT FORMAT  (CFS)
C..
       DO 975 M=1,DAYS_YR

            DO 957 K=1,NUMBER
               N=NODE(K)
               QDRN(N,YR,M)=ND(N,M)
               QIRR(N,YR,M)=NID1(N,M)
               QSEEP(N,YR,M)=NID2(N,M)
               IF(.NOT. NDFLG(N))THEN
                  IF((QDRN(N,YR,M)+QIRR(N,YR,M)+QSEEP(N,YR,M)).NE.0)THEN
                     NDFLG(N)=.TRUE.
                  ENDIF
               ENDIF
 957        CONTINUE

 975     CONTINUE
      ENDDO
C     end of the year cycle
c---------------------add groundwater output
       WRITE(56,2989)
 2989  FORMAT('GW_per_island.dss')
       DO ISL = 1,MAXISL
           WRITE(56,2910) ISL
 2910      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=GW-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(56,2911)GWF(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(56,2912)
      ENDDO
      WRITE(56,2913)
 2911 FORMAT(F10.2)
 2912 FORMAT('END')
 2913 FORMAT('FINISH')
      CLOSE(56)
c---------------------end of groundwater output
c---------------------add drainage_without_runoff output
       WRITE(58,2914)
 2914  FORMAT('drn_wo_ro_island.dss')
       DO ISL = 1,MAXISL
           WRITE(58,2915) ISL
 2915      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=DRN-WO-RO-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(58,2911)DRN_Y(ISL,M,YR)*TAF2CFS
               ENDDO
           ENDDO
           WRITE(58,2912)
      ENDDO
      WRITE(58,2913)
      CLOSE(58)
c---------------------end of drainage_without_runoff output
c---------------------add diversion_without_seepage output
       WRITE(60,2916)
 2916  FORMAT('div_wo_spg_island.dss')
       DO ISL = 1,MAXISL
           WRITE(60,2917) ISL
 2917      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=DIV-WO-SPG-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(60,2911)DIV_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(60,2912)
      ENDDO
      WRITE(60,2913)
      CLOSE(60)
c---------------------end of diversion without seepage output
c---------------------add seepage output
       WRITE(62,2918)
 2918  FORMAT('spg_island.dss')
       DO ISL = 1,MAXISL
           WRITE(62,2919) ISL
 2919      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=SPG-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(62,2911)SPG_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(62,2912)
      ENDDO
      WRITE(62,2913)
      CLOSE(62)
c---------------------end of seepage output
c---------------------add runoff output
       WRITE(64,2920)
 2920  FORMAT('RO_island.dss')
       DO ISL = 1,MAXISL
           WRITE(64,2921) ISL
 2921      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=RO-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(64,2911)RO_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(64,2912)
      ENDDO
      WRITE(64,2913)
      CLOSE(64)
c---------------------end of runoff output
c---------------------add leach drainage
       WRITE(66,2922)
 2922  FORMAT('LRD_island.dss')
       DO ISL = 1,MAXISL
           WRITE(66,2923) ISL
 2923      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=LRD-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(66,2911)LRD_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(66,2912)
      ENDDO
      WRITE(66,2913)
      CLOSE(66)
c---------------------end of leach drainage
c---------------------add leach applied water
       WRITE(68,2924)
 2924  FORMAT('LRA_island.dss')
       DO ISL = 1,MAXISL
           WRITE(68,2925) ISL
 2925      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=LRA-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(68,2911)LRA_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(68,2912)
      ENDDO
      WRITE(68,2913)
      CLOSE(68)
c---------------------end of leach applied water
c---------------------add temporary drainage
       WRITE(70,2926)
 2926  FORMAT('tempD_island.dss')
       DO ISL = 1,MAXISL
           WRITE(70,2927) ISL
 2927      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=tempDRNs-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(70,2911)tempD_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(70,2912)
      ENDDO
      WRITE(70,2913)
      CLOSE(70)
c---------------------end of temporary drainage
c---------------------add groundwater amount 2
       WRITE(72,2928)
 2928  FORMAT('GWAmount2_island.dss')
       DO ISL = 1,MAXISL
           WRITE(72,2929) ISL
 2929      FORMAT('A=DICU-ISLAND  B=',I4,
     &       '  C=GWA2-FLOW  D=01JAN2015  E=1DAY  F=DWR-BDO'/
     &       'CFS'/
     &       'PER-AVER'/
     &       '01OCT2015 2400')
           DO YR = iyr1,iyr2
               IF(MOD(YR,4).EQ.0) THEN
                   DAYS_YR = 366
               ELSE
                   DAYS_YR = 365
               ENDIF
               DO M=1,DAYS_YR
                   WRITE(72,2911)GWA2_Y(ISL,M,YR)* TAF2CFS
               ENDDO
           ENDDO
           WRITE(72,2912)
      ENDDO
      WRITE(72,2913)
      CLOSE(72)

      CLOSE(80)
      CLOSE(81)
      CLOSE(82)

C..
C  FORMAT STATEMENTS
C..
   10 FORMAT(//9X,366F8.0)
   11 FORMAT(I4,4X,366F8.0)
   12 FORMAT(I4,4X,F8.2)
   13 FORMAT(5X,366F8.0)
   15 FORMAT(8X,12F8.0)
   19 FORMAT(///)
   20 FORMAT(I4,1X,I4,2X,F10.2)
   21 FORMAT(1X,I4,1X,I4,2X,F10.2)
   22 FORMAT(I4,4X,2F10.2)
   23 FORMAT(2I4,366F8.2)
  996 FORMAT( 'DELTA ISLAND PRECIP AND TOTAL CONSUMPTIVE USE'/
     .'ECHO PRINT OF INPUT DATA FOR WATER YEAR ',I4)
  997 FORMAT('DWR DELTA MODEL NODAL CHANNEL DIVERSIONS BY ISLAND (CFS)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' ISL  NODE',
     .T8/)
  998 FORMAT('DWR DELTA MODEL NODAL DRAINAGE RETURNS BY ISLAND (CFS)
     &  WARNING: RUNOFF NOT INCLUDED'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' ISL  NODE',
     .T8/)
  999 FORMAT( 'DELTA ISLAND CHANNEL DIVERSIONS AND DRAINAGES'/
     .'ECHO PRINT OF INPUT DATA FOR WATER YEAR ',I4)
 1000 FORMAT(I4,4X,366F8.0,4X,F8.0)
 1001 FORMAT('DELTA ISLAND CHANNEL DIVERSIONS (ACRE-FEET)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' ISLAND',
     .T7,12(5X,A3),4X,'ISL TOTAL'/)
 1002 FORMAT(//'DELTA ISLAND DRAINAGE RETURNS (ACRE-FEET)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' ISLAND',
     .T7,12(5X,A3),4X,'ISL TOTAL'/)
 1003 FORMAT('DWRDSM DELTA MODEL NODAL CHANNEL DIVERSIONS (CFS)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' NODE',
     .T7,12(5X,A3),4X/)
 1004 FORMAT(//'DWRDSM DELTA MODEL NODAL DRAINAGE RETURNS (CFS)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' NODE',
     .T7,12(5X,A3)/)
 1005 FORMAT('TOTAL: MONTHLY AND ANNUAL DELTA DIVERSIONS (ACRE-FT)'/
     .       8X,366F8.0,2X,F10.0)
 1006 FORMAT('TOTAL: MONTHLY AND ANNUAL DELTA DRAINAGE (ACRE-FT)'/
     .       8X,366F8.0,2X,F10.0)
CNM---1007 FORMAT('DWR/RMA DELTA MODEL HYDROLOGY ENTRIES BY NODE (CFS)'/
CNM---.'DRAINAGES (INPUT) AND DIVERSIONS (OUTPUT) FOR ',A3,
CNM---.' - WATER YEAR ',I4 //' NODE      DRN	     DIV'/)
CNM---1008 FORMAT(I5,2F10.2)

 1007 FORMAT('DWRDSM DELTA MODEL HYDROLOGY ENTRIES BY NODE (CFS)'/,
     .'DRAINAGES (INPUT) AND DIVERSIONS (OUTPUT) FOR ',
     &' - WATER YEAR ',I4 /,
     .'CHANNEL DIVERSIONS ARE SEPARATED INTO IRRIGATION
     & DIVERSIONS AND SEEPAGE'/,
     .' NODE      DRN	   CH DIV  IRRIG DIV   SEEPAGE  '/)

 1008  FORMAT(I5,4F10.2)

 1009 FORMAT(I4,4X,366F8.2)
 1010 FORMAT(/'MONTHLY DELTA DIVERSIONS OVER ALL NODES (CFS)'/
     .       8X,366F8.0)
 1011 FORMAT(/'MONTHLY DELTA DRAINAGE RETURNS OVER ALL NODES (CFS)'/
     .       8X,366F8.0)
 1012 FORMAT(//'MONTHLY AND ANNUAL NET DELTA CD (ACRE-FT)'/
     .       8X,366F8.0,2X,F10.0)
 1013 FORMAT('DWRDSM DELTA MODEL NODAL CHANNEL DIVERSIONS (ACRE-FT)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' NODE',
     .T7,12(5X,A3)/)
 1014 FORMAT(//'DWRDSM DELTA MODEL NODAL DRAINAGE RETURNS (ACRE-FT)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' NODE',
     .T7,12(5X,A3)/)
 1015 FORMAT(//'MONTHLY NET DELTA CD OVER ALL NODES (AF)'/
     .       8X,366F8.0)
 1016 FORMAT(//'MONTHLY NET DELTA CD OVER ALL NODES (CFS)'/
     .       8X,366F8.0)
 1017 FORMAT('DWR DELTA MODEL NODAL DRAINAGE TDS CONCENTRATION (MG/L)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' NODE',
     .T7,12(5X,A3),4X/)
 1018 FORMAT('DWR DELTA MODEL NODAL DRAINAGE TDS CONCENTRATION (MG/L)'/
     .'AND DRAINAGE RETURN FLOW (CFS) FOR  ',A3,' - WATER YEAR ',I4 //
     .' NODE    MG/L      CFS'/)
 2018 FORMAT('DWR DELTA MODEL NODAL DRAINAGE CL CONCENTRATION (MG/L)'/
     .'AND DRAINAGE RETURN FLOW (CFS) FOR  ',A3,' - WATER YEAR ',I4 //
     .' NODE    MG/L      CFS'/)
 1019 FORMAT(I4,4X,366F8.0)
 1049 FORMAT(I4,4X,12F8.0)
 1020 FORMAT('MONTHLY DELTA DIVERSIONS OVER ALL NODES (AF)'/
     .       8X,366F8.0)
 1021 FORMAT('MONTHLY DELTA DRAINAGE RETURNS OVER ALL NODES (AF)'/
     .       8X,366F8.0)
 1022 FORMAT(I5,F10.0,F10.2)
C-----1023 FORMAT('DWR/RMA DELTA MODEL HYDROLOGY ENTRIES BY ISLAND (CFS)'/
C-----.'DRAINAGES (INPUT) AND DIVERSIONS (OUTPUT) FOR ',A3,
C-----.' - WATER YEAR ',I4 //' ISLAND CH. DEP    AG-RET    P - ET'/)
 1023 FORMAT('DWRDSM DELTA MODEL HYDROLOGY ENTRIES BY ISLAND (CFS)'/
     .'DRAINAGES (INPUT) AND DIVERSIONS (OUTPUT) FOR ',A3,
     .' - WATER YEAR ',I4 //' ISLAND CH. DEP    AG-RET    H(FT.)'/)
C..

CNM-10/04/95
 1025 FORMAT(//'DELTA ISLAND CHANNEL DIVERSIONS: AW+LW, NO SEEPAGE
     & INCLUDED (ACRE-FEET)'/
     .'FOR OCTOBER THRU SEPTEMBER - WATER YEAR ',I4 //' ISLAND',
     .T7,12(5X,A3))

 1024 FORMAT(I4,4X,366F8.0)

CNM-10/04/95

      STOP 'SUCCESS'
      END

