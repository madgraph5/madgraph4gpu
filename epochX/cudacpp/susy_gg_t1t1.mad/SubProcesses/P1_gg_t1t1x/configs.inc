C     Diagram 1
      DATA MAPCONFIG(1)/1/
      DATA TSTRATEGY(1)/2/
      DATA (IFOREST(I,-1,1),I=1,2)/4,3/
      DATA (SPROP(I,-1,1),I=1,1)/21/
      DATA TPRID(-1,1)/0/
C     Diagram 2
      DATA MAPCONFIG(2)/2/
      DATA TSTRATEGY(2)/2/
      DATA (IFOREST(I,-1,2),I=1,2)/1,3/
      DATA TPRID(-1,2)/1000006/
      DATA (SPROP(I,-1,2),I=1,1)/0/
      DATA (IFOREST(I,-2,2),I=1,2)/-1,4/
      DATA TPRID(-2,2)/21/
      DATA (SPROP(I,-2,2),I=1,1)/0/
C     Diagram 3
      DATA MAPCONFIG(3)/3/
      DATA TSTRATEGY(3)/2/
      DATA (IFOREST(I,-1,3),I=1,2)/1,3/
      DATA TPRID(-1,3)/2000006/
      DATA (SPROP(I,-1,3),I=1,1)/0/
      DATA (IFOREST(I,-2,3),I=1,2)/-1,4/
      DATA TPRID(-2,3)/21/
      DATA (SPROP(I,-2,3),I=1,1)/0/
C     Diagram 4
      DATA MAPCONFIG(4)/4/
      DATA TSTRATEGY(4)/2/
      DATA (IFOREST(I,-1,4),I=1,2)/1,4/
      DATA TPRID(-1,4)/1000006/
      DATA (SPROP(I,-1,4),I=1,1)/0/
      DATA (IFOREST(I,-2,4),I=1,2)/-1,3/
      DATA TPRID(-2,4)/21/
      DATA (SPROP(I,-2,4),I=1,1)/0/
C     Diagram 5
      DATA MAPCONFIG(5)/5/
      DATA TSTRATEGY(5)/2/
      DATA (IFOREST(I,-1,5),I=1,2)/1,4/
      DATA TPRID(-1,5)/2000006/
      DATA (SPROP(I,-1,5),I=1,1)/0/
      DATA (IFOREST(I,-2,5),I=1,2)/-1,3/
      DATA TPRID(-2,5)/21/
      DATA (SPROP(I,-2,5),I=1,1)/0/
C     Number of configs
      DATA MAPCONFIG(0)/5/
C     used fake id
      DATA FAKE_ID/7/
