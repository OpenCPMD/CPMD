echo "NAME='char_r1'     TYPE='character' KIND='(*)'       DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='len(data)*n'  MPI_TYPE='MPI_CHARACTER' SIZE_IN_BYTES='mp_char_in_bytes'"

echo "NAME='int4_r1'     TYPE='integer'   KIND='(int_4)'   DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_INTEGER' SIZE_IN_BYTES='mp_int4_in_bytes'"
echo "NAME='int4_r2'     TYPE='integer'   KIND='(int_4)'   DIMENSION='dimension(1,*)'       COPY='(1,1:n)'       ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_INTEGER' SIZE_IN_BYTES='mp_int4_in_bytes'"

echo "NAME='int8_r1'     TYPE='integer'   KIND='(int_8)'   DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_INTEGER8' SIZE_IN_BYTES='mp_int8_in_bytes'"

echo "NAME='real8_r1'    TYPE='real'      KIND='(real_8)'  DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_PRECISION'  SIZE_IN_BYTES='mp_double_in_bytes'"
echo "NAME='real8_r2'    TYPE='real'      KIND='(real_8)'  DIMENSION='dimension(1,*)'       COPY='(1,1:n)'       ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_PRECISION'  SIZE_IN_BYTES='mp_double_in_bytes'"
echo "NAME='real8_r3'    TYPE='real'      KIND='(real_8)'  DIMENSION='dimension(1,1,*)'     COPY='(1,1,1:n)'     ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_PRECISION'  SIZE_IN_BYTES='mp_double_in_bytes'"
echo "NAME='real8_r4'    TYPE='real'      KIND='(real_8)'  DIMENSION='dimension(1,1,1,*)'   COPY='(1,1,1,1:n)'   ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_PRECISION'  SIZE_IN_BYTES='mp_double_in_bytes'"
echo "NAME='real8_r5'    TYPE='real'      KIND='(real_8)'  DIMENSION='dimension(1,1,1,1,*)' COPY='(1,1,1,1,1:n)' ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_PRECISION' SIZE_IN_BYTES='mp_double_in_bytes'"

echo "NAME='complex8_r1' TYPE='complex'   KIND='(real_8)'  DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_COMPLEX' SIZE_IN_BYTES='mp_double_complex_in_bytes'"
echo "NAME='complex8_r2' TYPE='complex'   KIND='(real_8)'  DIMENSION='dimension(1,*)'       COPY='(1,1:n)'       ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_COMPLEX' SIZE_IN_BYTES='mp_double_complex_in_bytes'"
echo "NAME='complex8_r3' TYPE='complex'   KIND='(real_8)'  DIMENSION='dimension(1,1,*)'     COPY='(1,1,1:n)'     ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_COMPLEX' SIZE_IN_BYTES='mp_double_complex_in_bytes'"
echo "NAME='complex8_r4' TYPE='complex'   KIND='(real_8)'  DIMENSION='dimension(1,1,1,*)'   COPY='(1,1,1,1:n)'   ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_DOUBLE_COMPLEX' SIZE_IN_BYTES='mp_double_complex_in_bytes'"

echo "NAME='logical_r1'  TYPE='logical'   KIND=''          DIMENSION='dimension(*)'         COPY='(1:n)'         ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_LOGICAL' SIZE_IN_BYTES='mp_logical_in_bytes'"
echo "NAME='logical_r2'  TYPE='logical'   KIND=''          DIMENSION='dimension(1,*)'       COPY='(1,1:n)'       ARG_LEN_NAME='n'  MPI_LEN='n'  MPI_TYPE='MPI_LOGICAL' SIZE_IN_BYTES='mp_logical_in_bytes'"
