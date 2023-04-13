__global__ void DataPassGB(int *VtoC,int *CtoV,int *Receivedword, int *Interleaver,int Dev_ColumnDegree,int N,int NbBranch, int inter);
__global__ void CheckPassGB(int *DCtoV, int *VtoC,int M,int NbBranch,int Dev_RowDegree);
__global__ void APP_GB(int *Decide,int *CtoV,int *Receivedword,int *Interleaver,int Dev_ColumnDegree,int N,int M,int NbBranch);
__global__ void ComputeSyndrome(int *Decide,int *Mat,int Dev_RowDegree,int M, int *Dev_Syndrome);
