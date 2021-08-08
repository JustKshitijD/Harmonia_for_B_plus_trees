# include <bits/stdc++.h>
# define MAX 1000000000
using namespace std;

__device__ int k1_gpu=364; __device__ int k2_gpu=121;   // k1= 3^0+...+3^5 = 364= sizeof(nodes)
                          //  k2=364-3^5= sizeof(prefix_sum)
__device__ int m_gpu=2;                  // 1 key + 4 data items

float avg_tree_size=0.0;
int T_count=0;

__device__ int GS=3; __device__ int no_of_queries_in_warp=6;
int GS_cpu=3; int no_of_queries_in_warp_cpu=6;
int fanout=4; __device__ int fanout_gpu=8;             // so, each node can have max fanout-1 keys and fanout children
int nk=3;
__device__ int nk_gpu=3;


__device__ typedef struct Node
{
  int* key;      // keys are sorted
  int** data;
}node_gpu;

__device__ node_gpu** nodes_gpu;
__device__ int* prefix_sum_gpu;

int k1=364; int k2=121;

int m=2;
int sum_of_den=0;

bool flag=false;
int incremented_index=-1;


typedef struct Node_gpu
{
  int* key;      // keys are sorted
  int** data;
}node;


node** nodes;
int* prefix_sum;

int lb(int* arr, int N, int X)                     // binary search to find element, and give found index or index of next-highest element
{
	int mid;

	// Initialise starting index and
	// ending index
	int low = 0;
	int high = N;

	// Till low is less than high
	while (low < high) {
		mid = low + (high - low) / 2;

		// If X is less than or equal
		// to arr[mid], then find in
		// left subarray
		if (X <= arr[mid]) {
			high = mid;
		}

		// If X is greater arr[mid]
		// then find in right subarray
		else {
			low = mid + 1;
		}
	}

	// Return the lb index
	return low;
}


void init(node* nn)
{
  nn->key=(int*)malloc(sizeof(int)*nk);                // init each node with nk keys
  nn->data=(int**)malloc(sizeof(int*)*nk);             // and nk data arrays, each having m integers
  for(int i=0;i<nk;i++)
  {
    nn->data[i]=(int*)malloc(sizeof(int)*m);
  }

  for(int i=0;i<nk;i++)
    nn->key[i]=MAX;

    for(int i=0;i<nk;i++)
    {
      for(int j=0;j<m;j++)
      {
        nn->data[i][j]=MAX;
      }
    }
}

__device__ void init_gpu(node_gpu* nn)
{
    nn->key=(int*)malloc(sizeof(int)*nk_gpu);                // init each node with 3 keys
  nn->data=(int**)malloc(sizeof(int*)*nk_gpu);             // and 3 data arrays, each having m integers
  for(int i=0;i<nk_gpu;i++)
  {
    nn->data[i]=(int*)malloc(sizeof(int)*m_gpu);
  }

  for(int i=0;i<nk_gpu;i++)
    nn->key[i]=MAX;

    for(int i=0;i<nk_gpu;i++)
    {
      for(int j=0;j<m_gpu;j++)
      {
        nn->data[i][j]=MAX;
      }
    }
}

__global__ void iii()                                 // initializing nodes array in GPU
{
  nodes_gpu=(node_gpu**)malloc(sizeof(node_gpu*)*k1_gpu);
  for(int i=0;i<k1_gpu;i++)
    nodes_gpu[i]=NULL;
}

__global__ void cpy_key(int* key,int *i2)             // copy key to a node in GPU
{
    int i=*i2;
    printf("In cpy_key, i: %d\n",i);
    nodes_gpu[i]=(node_gpu*)malloc(sizeof(node_gpu));
    init_gpu(nodes_gpu[i]);

    for(int h=0;h<nk_gpu;h++)
    {
        nodes_gpu[i]->key[h]=key[h];
    }

}

__global__ void cpy_data(int* data,int* i2, int* ii2)    // ith node in nodes[]; iith key in that node
{
    int i=*i2; int ii=*ii2;
    nodes_gpu[i]->data[ii]=data;
    for(int h=0;h<m_gpu;h++)
    {
        (nodes_gpu[i]->data)[ii][h]=data[h];
    }
}

__global__ void print()
{
  printf("In print(), k1_gpu: %d\n",k1_gpu);
    for(int i=0;i<k1_gpu && nodes_gpu[i]!=NULL;i++)
    {
        printf("nodes_gpu[%d]->key: \n",i);
        for(int j=0;j<nk_gpu;j++)
        {
            printf("%d ",(nodes_gpu[i]->key)[j]);
        }
        printf("\n");

        printf("nodes_gpu[%d]->data: \n",i);
        for(int j=0;j<nk_gpu;j++)
        {
            for(int kk=0;kk<m_gpu;kk++)
            {
              printf("%d ",(nodes_gpu[i]->data)[j][kk]);
            }
            printf("\n");
        }
        printf("\n");
    }
}


int search(int key)                               // search for index of node in nodes[], having the given key
{
  int i=0; //node* ptr; ptr=nodes[0];
  int j=0;

  while(true)
  {
    //ptr=nodes[i];
    j=i;
    if(nodes[i]==NULL)                     // no node inserted yet
      return 0;
    int* key_arr=nodes[i]->key;
    int ind=lb(key_arr,nk,key);
    cout<<"i: "<<i<<"; ind: "<<ind<<" ; prefix_sum[i]: "<<prefix_sum[i]<<endl;
    if(prefix_sum[i]==MAX)                  // reached a leaf
      break;
    i=prefix_sum[i]+ind;
    if(key_arr[ind]==key)
      i++;
  }

  int* key_arr=nodes[j]->key;
  int ind=lb(key_arr,nk,key);
  //cout<<"SEARCH"<<endl;
  if(key_arr[ind]==key)
  {
    int** dd=nodes[j]->data;
    for(int r=0;r<m;r++)
        cout<<dd[ind][r]<<" ";
    cout<<endl;
  }

  return j;
}

void range_query(int k1,int k2)                   // scan leaves from first key to last key
{
  int j=search(k1);
  int* key_arr=nodes[j]->key;
  int ind=lb(key_arr,nk,k1);
  //cout<<"SEARCH"<<endl;
  cout<<"range_query: "<<endl;
  if(key_arr[ind]==k1)
  {
    cout<<"In"<<endl;
    ind++;
    int** dd=nodes[j]->data;

    for(;ind<nk;ind++)
    {
      if(key_arr[ind]==MAX)
        break;
      for(int r=0;r<m;r++)
          cout<<dd[ind][r]<<" ";
     }
    cout<<endl;
  }

  j++; bool f2=false;

  while(nodes[j]!=NULL && f2==false)
  {
    int* key_arr=nodes[j]->key;
    int** dd=nodes[j]->data;
    int ind=0;
    for(;ind<nk ;ind++)
    {
      if(key_arr[ind]==MAX)
        break;
      if(dd[ind][0]>k2)
      {
        f2=true;
        break;
      }
      for(int r=0;r<m;r++)
          cout<<dd[ind][r]<<" ";
      cout<<endl;
    }
    j++;
  }
}

void put_in_middle(node* nn,int pos, int key, int* dd)   // put key and its data dd at position pos in nn
{
  int arr[nk]; int dt[nk][m];
  for(int i=0;i<nk;i++)
    arr[i]=nn->key[i];

  for(int i=0;i<nk;i++)
  {
    for(int j=0;j<m;j++)
    {
      dt[i][j]=nn->data[i][j];
    }
  }

  nn->key[pos]=key;
  for(int i=pos+1;i<nk;i++)
    nn->key[i]=arr[i-1];

  for(int j=0;j<m;j++)
    nn->data[pos][j]=dd[j];

    for(int i=pos+1;i<nk;i++)
    {
      for(int j=0;j<m;j++)
      {
        nn->data[i][j]=dt[i-1][j];
      }
    }
}

void put_in_nodes(int ind,node* nn2)            // put nn2 in nodes at ind position
{
  node* g=nodes[ind]; nodes[ind]=nn2; int i=ind+1;
  for(i=ind+1;i<k1 && nodes[i]!=NULL;i++)
  {
    node* h=nodes[i];
    nodes[i]=g;
    g=h;
  }
  nodes[i]=g;
}

void put_in_array(int* arr, int sz, int x, int ind)      // put value x at index ind in array arr. sz is size of array
{
  int g=arr[ind]; arr[ind]=x; int i=ind+1;
  for(i=ind+1;i<k1 && arr[i]!=MAX;i++)
  {
    int h=arr[i];
    arr[i]=g;
    g=h;
  }
  arr[i]=g;
}

int parent(int ind)                           // give index of parent of node ind
{
  if(prefix_sum[0]==MAX || ind==0)                      // no element in prefix_sum
    return -1;

  int p=lb(prefix_sum,k2,ind);
  //cout<<"pp: "<<p<<"; prefix_sum[p]: "<<prefix_sum[p]<<endl;
  if(prefix_sum[p]!=ind)
    p--;
  return p;
}


int insert_internal(node* nn,int ind,int key,int* dd,int orig)    // insert key in node* nn
{
  cout<<"start ind: "<<ind<<"; start key: "<<key<<endl;
  cout<<"Keys of start node nn: "<<endl;
  for(int i=0;i<nk;i++)
    cout<<nn->key[i]<<" ";
  cout<<endl;

  int pos=lb(nn->key,nk,key);
  if(nn->key[nk-1]==MAX)                  // space there in node
  {
      cout<<"Inserting normally"<<endl;
      put_in_middle(nn,pos,key,dd);
      return 0;

  }
  else                                    // node is full
  {
    cout<<"Node is full"<<endl;
    for(int i=0;i<nk;i++)
    {
      cout<<nn->key[i]<<"# "<<endl;
    }
    cout<<"key: "<<key<<endl;

    int pos=lb(nn->key,nk,key); node* nn2; int* dd_mid;
    dd_mid=(int*)malloc(m*sizeof(int));
   int middle_element;  // middle_element goes to top next time
   cout<<"pos: "<<pos<<endl;

    if(prefix_sum[ind]==MAX)               // leaf
    {
      int mid=(nk+1)/2;
      cout<<"mid: "<<mid<<endl;

      if(pos<mid)
      {
        nn2=(node*)malloc(sizeof(node));
        init(nn2); middle_element=nn->key[mid-1]; // mid-1 goes to top
        for(int j=0;j<m;j++)
          dd_mid[j]=nn->data[mid-1][j];

        for(int i=mid-1;i<nk;i++)
        {
          nn2->key[i-(mid-1)]=nn->key[i];
        }

        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-(mid-1)][j]=nn->data[i][j];
          }
        }

        for(int i=mid-1;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }

        put_in_middle(nn,pos,key,dd);

      }
      else
      {
        if(pos==mid)
        {
          middle_element=key;
          for(int j=0;j<m;j++)
            dd_mid[j]=dd[j];
        }
        else
        {
          middle_element=nn->key[mid];
          cout<<"middle_element: "<<middle_element<<endl;
          for(int j=0;j<m;j++)
          {
            cout<<"j: "<<j<<endl;
            dd_mid[j]=nn->data[mid][j];
          }
        }

        nn2=(node*)malloc(sizeof(node));
        init(nn2);

        for(int i=mid;i<nk;i++)
        {
          nn2->key[i-mid]=nn->key[i];
        }
        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-mid][j]=nn->data[i][j];
          }
        }
        put_in_middle(nn2,pos-mid,key,dd);

        cout<<"mid2: "<<mid<<endl;
        for(int i=mid;i<nk;i++)
        {
          cout<<"i: "<<i<<endl;
          nn->key[i]=MAX;
        }

        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }
      }

    }
    else
    {
      int mid=(nk+1)/2;

      if(pos<mid)
      {
        nn2=(node*)malloc(sizeof(node));
        init(nn2); middle_element=nn->key[mid-1]; // mid-1 goes to top
        for(int j=0;j<m;j++)
          dd_mid[j]=nn->data[mid-1][j];

        for(int i=mid;i<nk;i++)                  // leaving mid-1
        {
          nn2->key[i-(mid)]=nn->key[i];
        }


        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn2->data[i-(mid)][j]=nn->data[i][j];
          }
        }

        for(int i=mid-1;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid-1;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }

        put_in_middle(nn,pos,key,dd);

      }
      else
      {
        if(pos==mid)
        {
          middle_element=key;
          for(int j=0;j<m;j++)
            dd_mid[j]=dd[j];
          nn2=(node*)malloc(sizeof(node));
          init(nn2);

          for(int i=mid;i<nk;i++)
          {
            nn2->key[i-mid]=nn->key[i];
          }
          for(int i=mid;i<nk;i++)
          {
            for(int j=0;j<m;j++)
            {
              nn2->data[i-mid][j]=nn->data[i][j];
            }
          }

        }
        else
        {
          middle_element=nn->key[mid];
          for(int j=0;j<m;j++)
            dd_mid[j]=nn->data[mid][j];
          nn2=(node*)malloc(sizeof(node));
          init(nn2);

          for(int i=mid+1;i<nk;i++)
          {
            nn2->key[i-(mid+1)]=nn->key[i];
          }
          for(int i=mid+1;i<nk;i++)
          {
            for(int j=0;j<m;j++)
            {
              nn2->data[i-(mid+1)][j]=nn->data[i][j];
            }
          }
          put_in_middle(nn2,pos-(mid+1),key,dd);
        }

        for(int i=mid;i<nk;i++)
        {
          nn->key[i]=MAX;
        }
        for(int i=mid;i<nk;i++)
        {
          for(int j=0;j<m;j++)
          {
            nn->data[i][j]=MAX;
          }
        }
      }
    }

    put_in_nodes(ind+1,nn2);

    if(prefix_sum[ind]==MAX)               // leaf
    {
      int xx=0; int p=parent(ind);
      if(p==-1 || nodes[p]==NULL)     // leaf and parent is null means first node to be filled, and now, prefix_sum will get its first entry
      {
        cout<<endl;
        cout<<"Leaf with parent NULL:- ind: "<<ind<<"; p: "<<p<<endl;
        prefix_sum[0]=1;
        node* nn3; nn3=(node*)malloc(sizeof(node));
        init(nn3);
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_middle(nn3,0,middle_element,dd_mid);
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_nodes(0,nn3);
        return 0;
      }
      else
        xx=insert_internal(nodes[p],p,middle_element,dd_mid,ind);     // dd is NULL means we won't insert anything

      ind+=xx;
      cout<<"Leaf with parent non-null:- ind+=xx: "<<ind<<"; p: "<<p<<endl;


      cout<<endl;
      p=parent(ind);
      int ind2=prefix_sum[p];
      cout<<"prefix_sum: "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      for(int i=0;i<k1;i++)
      {
        if(nodes[i]!=NULL)
        {
          cout<<"nodes["<<i<<"]->key: "<<endl;
          for(int j=0;j<nk;j++)
          {
            cout<<(nodes[i]->key)[j]<<" ";
          }
          cout<<endl;
          cout<<"nodes["<<i<<"]->data: "<<endl;
          for(int j=0;j<nk;j++)
          {
            cout<<"data["<<j<<"]- "<<endl;
            for(int f=0;f<m;f++)
            {
              cout<<(nodes[i]->data)[j][f]<<" ";
            }
            cout<<endl;
          }
          cout<<endl;
        }
      }

      cout<<"p: "<<p<<" ; ind: "<<ind<<"; ind2: "<<ind2<<endl;

      cout<<"prefix_sum in full leaf whose parent is not null, before: "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      //if(flag==false)
      //{
      cout<<"incremented_index: "<<incremented_index<<endl;
        for(int i=0;i<k2;i++)
        {
          if(prefix_sum[i]==MAX)
            break;

          if(prefix_sum[i]>ind2)
          {
            cout<<"flag: "<<flag<<endl;
            if(!(prefix_sum[i]==incremented_index && flag==true))
              prefix_sum[i]+=1;
          }
        }
      //}

      cout<<"prefix_sum in full leaf whose parent is not null,after : "<<endl;
      for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
      {
        cout<<prefix_sum[i]<<" ";
      }
      cout<<endl;

      return xx+1;
    }
    else                                                    // non-leaf
    {
      int x=0; int no_of_keys_of_v1=0; int xx=0;
      for(int i=0;i<nk;i++)
      {
        if(nn->key[i]!=MAX)
          no_of_keys_of_v1++;
        else
            break;
      }

      int p=parent(ind);
      if(p==-1 || nodes[p]==NULL)                              // parent is NULL
      {
        cout<<endl;
        cout<<"In non-leaf with parent NULL, ind: "<<ind<<"; parent(ind): "<<p<<endl;
        put_in_array(prefix_sum,k2,1,0);
        x=2;

        node* nn3; nn3=(node*)malloc(sizeof(node));
        init(nn3);
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_middle(nn3,0,middle_element,dd_mid);
        //cout<<"middle_element: "<<middle_element<<endl;
        put_in_nodes(0,nn3);

        int pf_2=-1;
        int ss=prefix_sum[1]+2;
        cout<<"ss: "<<ss<<" ; orig+2: "<<orig+2<<"; ind: "<<ind<<endl;
        if(ss>orig+2)
        {
          pf_2=ss-2+no_of_keys_of_v1;
          flag=false;
        }
        else
        {
          pf_2=ss-2+no_of_keys_of_v1+1;
          flag=true;
        }

        cout<<"pf_2: "<<pf_2<<endl;
        incremented_index=pf_2+2;;
        cout<<"incremented_index in orig+2: "<<incremented_index<<endl;
        put_in_array(prefix_sum,k2,pf_2,2);

        cout<<"prefix_sum before "<<endl;
        for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
        {
          cout<<prefix_sum[i]<<" ";
        }
        cout<<endl;

        for(int i=1;i<k2;i++)
        {
          if(prefix_sum[i]==MAX)
            break;

          if(prefix_sum[i]>ind)
            prefix_sum[i]+=x;
        }

        cout<<"prefix_sum after: "<<endl;
        for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
        {
          cout<<prefix_sum[i]<<" ";
        }
        cout<<endl;

        return 2;
      }
      else
      {
        int p=parent(ind);
        int xx=insert_internal(nodes[p],p,middle_element,dd_mid,ind);
        ind+=xx;
        x=1;
      }
      cout<<endl;
      cout<<"In non-leaf with parent non-null, ind+=xx: "<<ind<<"; parent(ind): "<<p<<endl;

      for(int i=0;i<k2;i++)
      {
        if(prefix_sum[i]==MAX)
          break;

        if(prefix_sum[i]>ind)
          prefix_sum[i]+=x;
      }

      int pf_2=-1;
      int ss=prefix_sum[ind];
      cout<<"ss: "<<ss<<" ; orig+xx+1: "<<orig+xx+1<<"; ind: "<<ind<<"; xx: "<<xx<<"; orig: "<<orig<<endl;
      if(ss>orig+xx+1)
      {
        pf_2=ss+no_of_keys_of_v1;
        flag=false;
      }
      else
      {
        pf_2=ss+no_of_keys_of_v1+1;
        flag=true;
      }

      cout<<"pf_2: "<<pf_2<<endl;
      incremented_index=pf_2;
      cout<<"incremented_index in orig+ind: "<<incremented_index<<endl;
      put_in_array(prefix_sum,k2,pf_2,ind+1);

      return (xx+1);
    }

  }
}

void insert(int* dd)               // insert key and data dd into nodes and update prefix_sum
{
  int key=dd[0];
  int ind=search(key);             // gives index of leaf in nodes[], where key can be inserted
  node* nn=nodes[ind];

  if(nn==NULL)                      // nodes[0] is NULL
  {
    nn=(node*)malloc(sizeof(node));
    init(nn);
    nn->key[0]=key;
    for(int j=0;j<m;j++)
    {
      nn->data[0][j]=dd[j];
    }

    nodes[0]=nn;

  }
  else                            // we have the leaf node nn to which original function can be applied
  {
    insert_internal(nn,ind,key,dd,-1);
  }

  cout<<"prefix_sum: "<<endl;
  for(int i=0;i<k2 && prefix_sum[i]!=MAX;i++)
  {
    cout<<prefix_sum[i]<<" ";
  }
  cout<<endl;

  cout<<"After inserting key "<<dd[0]<<": "<<endl;
  for(int i=0;i<k1;i++)
  {
    if(nodes[i]!=NULL)
    {
      cout<<"nodes["<<i<<"]->key: "<<endl;
      for(int j=0;j<nk;j++)
      {
        cout<<(nodes[i]->key)[j]<<" ";
      }
      cout<<endl;
      cout<<"nodes["<<i<<"]->data: "<<endl;
      for(int j=0;j<nk;j++)
      {
        cout<<"data["<<j<<"]- "<<endl;
        for(int f=0;f<m;f++)
        {
          cout<<(nodes[i]->data)[j][f]<<" ";
        }
        cout<<endl;
      }
      cout<<endl;
    }
  }

}


int getMax(int* arr, int n)      // get max value in array
{
  cout<<"In get_max"<<endl;
  cout<<"n: "<<n<<endl;
  for(int i=0;i<n;i++)
    cout<<arr[i]<<" ";
  cout<<endl;
	int mx = arr[0];
	for (int i = 1; i < n; i++)
		if (arr[i] > mx)
			mx = arr[i];
	return mx;
}

// A function to do counting sort of arr[] according to
// the digit represented by exp.
void countSort(int arr[], int n, int exp)
{
	int output[n]; // output array
	int i, count[10] = { 0 };

	// Store count of occurrences in count[]
	for (i = 0; i < n; i++)
		count[(arr[i] / exp) % 10]++;

	// Change count[i] so that count[i] now contains actual
	// position of this digit in output[]
	for (i = 1; i < 10; i++)
		count[i] += count[i - 1];

	// Build the output array
	for (i = n - 1; i >= 0; i--) {
		output[count[(arr[i] / exp) % 10] - 1] = arr[i];
		count[(arr[i] / exp) % 10]--;
	}

	// Copy the output array to arr[], so that arr[] now
	// contains sorted numbers according to current digit
	for (i = 0; i < n; i++)
		arr[i] = output[i];
}

void radixsort(int* arr, int d, int N)          // only sort based on Nth bit from LSB to MSB
{
  cout<<"d: "<<d<<endl;
	// Find the maximum number to know number of digits
	int m = getMax(arr, d);
  printf("m: %d; N:%d\n",m,N );
  //cout<<"m: "<<m<<" ; N: "<<N<<endl;

	for (int exp = N; m / exp > 0; exp *= 10)
		countSort(arr, d, exp);
}

void pre_sort(int* c, int d)
{
	int B=7; int K=16;
  cout<<"In pre_sort"<<endl;
  int T=0;
  for(T=0;T<k1;T++)
  {
    if(nodes[T]==NULL)
      break;
  }
 avg_tree_size+=T;     // T is number of nodes present
 T_count++;

  cout<<"size of tree T: "<<T<<endl;
	int N=B-floor(log2(2^B/(T)*K));
  cout<<"N: "<<N<<endl;
	radixsort(c,d,N);
}

__global__ void srch(int* c, int* ll, int* query_index_arr,int* max_comp, int* d, int* qq)  // GS=3, no_of_queries_in_warp=6
{
	int id=blockIdx.x*blockDim.x+threadIdx.x; int q=*qq;
  printf("*ll: %d\n",*ll);
  printf("nk_gpu: %d\n",nk_gpu);
  int ii=(id/32)*no_of_queries_in_warp+(id%32)/GS;
  if(id==0)
  {
    printf("c:- \n");
    for(int i=0;i<*d;i++)
    {
      printf("%d ",c[i]);
    }
    printf("\n");
  }

	if(id<*ll && id%32<no_of_queries_in_warp*GS && ii<*d)
	{
		int cmp=0;                      // 3 id's processing 1 query
		int key=c[ii];                                      // id of same query get same key
		//int xx;
		int i=0; //node* ptr; ptr=nodes[0];
	  int j=0;
    printf("id: %d; ii: %d\n",id,ii);
    if(id==3)
    {
      printf("prefix_sum_gpu: \n");
      for(int i=0;i<k2_gpu;i++)
      {
        printf("%d ",prefix_sum_gpu[i]);
      }
      printf("\n");
    }

	  while(true)
	  {
      for(int r=0;r<q;r++)
        query_index_arr[r]=-1;
	    //ptr=nodes[i];
      printf("In i: %d\n",i);
	    j=i;
	    if(nodes_gpu[i]==NULL)                     // no node inserted yet
	      return;
	    int* key_arr=nodes_gpu[i]->key;  int ke_id=-1;

        printf("key_arr for id=%d: \n",id);
        for(int r=0;r<nk_gpu;r++)
        {
          printf("key_arr[%d]- %d ",r,key_arr[r]);
        }
        printf("\n");


			if(key_arr[0]>key)
			{
        printf("1st- id: %d, key: %d, ii:%d\n",id,key,ii);
				ke_id=0; query_index_arr[ii]=0;
        cmp+=1;
			}
			else if(key_arr[nk_gpu-1]<key)
			{
        printf("2nd- id: %d, key: %d, ii:%d\n",id,key,ii);
				ke_id=nk_gpu; query_index_arr[ii]=nk_gpu;
        cmp+=1;
			}
			else
			{
        printf("3rd- id: %d, key: %d, ii:%d\n",id,key,ii);

				for(int k=id%GS;k<nk_gpu && query_index_arr[ii]==-1 ;k+=GS)
				{
          printf("id: %d, k:%d, ii:%d\n",id,k,ii);
          cmp+=GS;
					if(key_arr[k]==key)
					{
						ke_id=k;
						query_index_arr[ii]=ke_id;
            //if(prefix_sum_gpu[prefix_sum_gpu[i]+k]==MAX)
            if(prefix_sum_gpu[i]==MAX)
              printf("%d and %d are EQUAL!!!!!\n",key,key_arr[k]);
						break;
					}
					else if(k+1<nk_gpu && key_arr[k]<key && key_arr[k+1]>key)
					{
						ke_id=k+1;
						query_index_arr[ii]=ke_id;
            printf("In id:%d ;ke_id: %d\n",id,ke_id);
						break;
					}
          printf("id:%d ;ke_id: %d\n",id,ke_id);
				}
			}

      printf("query_index_arr[%d]= %d\n",ii,query_index_arr[ii]);
			int ind=query_index_arr[ii];       // ind,i same for all threads of query
	    //int ind=lb(key_arr,nk,key);

      printf("i: %d; ind: %d; prefix_sum_gpu[i]: %d\n",i,ind,prefix_sum_gpu[i]);
	    if(prefix_sum_gpu[i]==MAX)                  // reached a leaf
	      break;
	    i=prefix_sum_gpu[i]+ind;
	    if(key_arr[ind]==key)
	      i++;

      printf("Next going to i=%d\n",i);
      printf("\n");
                        			// threads that change a particular query_index_arr[ii] belong to same warp. So, barriers not needed anywhere.
	  }

	  int* key_arr=nodes_gpu[j]->key; query_index_arr[ii]=-1; int ke_id=-1;

		for(int k=id%GS; k<nk_gpu && query_index_arr[ii]==-1 ;k+=GS)
		{
			if(key_arr[k]==key)
			{
				ke_id=k;
				query_index_arr[ii]=ke_id;
				break;
			}
		}
	  //cout<<"SEARCH"<<endl;
	  if(key_arr[query_index_arr[ii]]==key)
	  {
	    int** dd=nodes_gpu[j]->data;
	    for(int r=0;r<m_gpu;r++)
	        printf("%d\n",dd[query_index_arr[ii]][r]);
	    printf("\n");
	  }

		atomicMax(max_comp,cmp);    // among all queries, max_comp is # of comparison steps of the querry that makes the warp do maximum comparisons
		printf("max_comp: %d\n",*max_comp);
	}
}

__global__ void init_pref_sum(int* pref)                 // Initialise prefix_sum_gpu using pref
{
  prefix_sum_gpu=(int*)malloc(sizeof(int)*k2_gpu);
  for(int i=0;i<k2_gpu;i++)
  {
    prefix_sum_gpu[i]=pref[i];
  }
}

void update()                                          // copies nodes[] and prefix_sum from CPU to GPU
{
  //int D=max_diff(b,q);
  int* pref;
  cudaMalloc(&pref,sizeof(int)*k2);
  cudaMemcpy(pref,prefix_sum,sizeof(int)*k2,cudaMemcpyHostToDevice);
  init_pref_sum<<<1,1>>>(pref);

  iii<<<1,1>>>();
  cudaDeviceSynchronize();

  printf("nodes_gpu initialized\n");

 for(int i=0;i<k1 && nodes[i]!=NULL;i++)
 {
     int* gpu_key; int* i2;  cudaMalloc(&gpu_key,sizeof(int)*nk); cudaMalloc(&i2,sizeof(int));
     cudaMemcpy(gpu_key,nodes[i]->key,sizeof(int)*nk,cudaMemcpyHostToDevice);
     cudaMemcpy(i2,&i,sizeof(int),cudaMemcpyHostToDevice);

     cpy_key<<<1,1>>>(gpu_key,i2);
     cudaDeviceSynchronize();
     cout<<"Key "<<i<<" copied"<<endl;

     for(int ii=0;ii<nk;ii++)
     {
       int* ii2;  cudaMalloc(&ii2,sizeof(int)); cudaMemcpy(ii2,&ii,sizeof(int),cudaMemcpyHostToDevice);

       //int* j2; cudaMalloc(&j2,sizeof(int)); cudaMemcpy(j2,&j,sizeof(int),cudaMemcpyHostToDevice);
       int* gpu_data; cudaMalloc(&gpu_data,sizeof(int)*m);
       cudaMemcpy(gpu_data,(nodes[i]->data)[ii],sizeof(int)*m,cudaMemcpyHostToDevice);

       cpy_data<<<1,1>>>(gpu_data,i2,ii2);
       cudaDeviceSynchronize();
     }

     cout<<"Data "<<i<<" copied"<<endl;

 }

  cout<<"Printing nodes_gpu on GPU"<<endl;
  print<<<1,1>>>();
  cudaDeviceSynchronize();
}


void func(int** a, int *b, int q)                                           // a has data of queries; b has index of inserts in a
{
  cout<<"b:-"<<endl;
  for(int i=0;i<q;i++)
  {
    cout<<b[i]<<" ";
  }
  cout<<endl; int vald=0;

	for(int i=0;i<q-1 && b[i]!=MAX && b[i+1]!=MAX;i++)         // 1 index of a has the element to be inserted. 0th has the number 2. b[i] is the index in a, where node to be inserted is present
	{
    vald++;                     // vald=len(b)-1
		int d=b[i+1]-b[i]-1;        // d queries
    cout<<"i: "<<i<<" ; q: "<<q<<endl;
    cout<<"b[i]: "<<b[i]<<" ; b[i+1]: "<<b[i+1]<<endl;
                                               // d queries done in paallel
    int* aa; aa=(int*)malloc(sizeof(int)*m);  cout<<"Inserting"<<endl;

    for(int j=1;j<=m;j++)
    {
      aa[j-1]=a[b[i]][j];
      cout<<aa[j-1]<<" ";
    }
    cout<<endl;

    cout<<"Inserting b["<<i<<"]"<<endl;
    flag=false;
		insert(aa);
    update();                           // call update everytime after insert to copy nodes to nodes_gpu


    cout<<"Back from insert"<<endl;
    int* c; c=(int*)malloc(sizeof(int)*d); int cnt=0;

    for(int j=b[i]+1;j<b[i+1];j++)                              // int* c=[a[b[i]] to a[b[i+1]]];
    {
      c[j-b[i]-1]=a[j][1]; cnt++;
    }
    cout<<"c: "<<endl;
    for(int j=0;j<d;j++)
    {
      cout<<c[j]<<" ";
    }
    cout<<endl;
    cout<<"cnt: "<<cnt<<endl;


    if(cnt>0)
    {
  		pre_sort(c,d);
      cout<<"c after pre-sorting: "<<endl;
      for(int j=0;j<d;j++)
      {
        cout<<c[j]<<" ";
      }
      cout<<endl;


  		int* query_index_arr; int* q_cpu;   int* max_comp;int* yy;   int* ll_cpu; int* ll_gpu;   int* c_gpu;  int* max_cmp_cpu;   int* d_gpu; int* qq;
  		q_cpu=(int*)malloc(sizeof(int)*d); yy=(int*)malloc(sizeof(int)); ll_cpu=(int*)malloc(sizeof(int)); max_cmp_cpu=(int*)malloc(sizeof(int));
  		for(int j=0;j<d;j++)
  			q_cpu[j]=-1;

  		*yy=0; *ll_cpu=ceil(d*1.0/no_of_queries_in_warp_cpu)*32;   // ll threads required for the d queries

  		cudaMalloc(&query_index_arr,sizeof(int)*d); cudaMalloc(&max_comp,sizeof(int));
  		cudaMemcpy(query_index_arr,q_cpu,sizeof(int)*d,cudaMemcpyHostToDevice);

  		cudaMalloc(&max_comp,sizeof(int));
  		cudaMemcpy(max_comp,yy,sizeof(int),cudaMemcpyHostToDevice);

      cout<<"*ll_cpu: "<<*ll_cpu<<endl;
  		cudaMalloc(&ll_gpu,sizeof(int));
  		cudaMemcpy(ll_gpu,ll_cpu,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&d_gpu,sizeof(int));
      cudaMemcpy(d_gpu,&d,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&qq,sizeof(int));
      cudaMemcpy(qq,&q,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&c_gpu,sizeof(int)*d);
  		cudaMemcpy(c_gpu,c,sizeof(int)*d,cudaMemcpyHostToDevice);

      cout<<"Before srch call"<<endl;
      cout<<"ceil(d*1.0/no_of_queries_in_warp_cpu): "<<ceil(d*1.0/no_of_queries_in_warp_cpu)<<endl;

  		srch<<<ceil(d*1.0/no_of_queries_in_warp_cpu),32>>>(c_gpu,ll_gpu,query_index_arr,max_comp,d_gpu,qq);   // max_comp is count of max_comparisons made in this search
  		cudaDeviceSynchronize();                                    // 1 block serves no_of_queries_in_warp queries
      cout<<"After srch"<<endl;
  		cudaMemcpy(max_cmp_cpu,max_comp,sizeof(int),cudaMemcpyDeviceToHost);

  		int product=(GS_cpu)*(*max_cmp_cpu);
      sum_of_den+=product;
      printf("max_cmp_cpu: %d\n",*max_cmp_cpu);
    }
	}

  int* aa; aa=(int*)malloc(sizeof(int)*m);
  cout<<"Inserting"<<endl;
  for(int j=1;j<=m;j++)
  {
    aa[j-1]=a[b[vald]][j];
    cout<<aa[j-1]<<" ";
  }
  cout<<endl;

  cout<<"Inserting b["<<vald<<"]"<<endl;
  flag=false;
  insert(aa);
  update();

  cout<<"vald: "<<vald<<endl; cout<<"b[vald]: "<<b[vald]<<"; b[vald+1]: "<<b[vald+1]<<endl;
  cout<<"q: "<<q<<endl;
  if(b[vald]<q)
  {
    int d=q-b[vald]-1;
    cout<<"d: "<<d<<endl;
    int* c; c=(int*)malloc(sizeof(int)*d); int cnt=0; int i=vald;

    for(int j=b[i]+1;j<q;j++)                              // int* c=[a[b[i]] to a[b[i+1]]];
    {
      c[j-b[i]-1]=a[j][1]; cnt++;
    }
    cout<<"c: "<<endl;
    for(int j=0;j<d;j++)
    {
      cout<<c[j]<<" ";
    }
    cout<<endl;
    cout<<"cnt: "<<cnt<<endl;

    if(cnt>0)
    {
      cout<<"Before pre_sorting"<<endl;
  		pre_sort(c,d);
      cout<<"c after pre-sorting: "<<endl;
      for(int j=0;j<d;j++)
      {
        cout<<c[j]<<" ";
      }
      cout<<endl;

  		int* query_index_arr; int* q_cpu;   int* max_comp;int* yy;   int* ll_cpu; int* ll_gpu;   int* c_gpu;  int* max_cmp_cpu;int* d_gpu;  int* qq;
  		q_cpu=(int*)malloc(sizeof(int)*d); yy=(int*)malloc(sizeof(int)); ll_cpu=(int*)malloc(sizeof(int)); max_cmp_cpu=(int*)malloc(sizeof(int));
  		for(int j=0;j<d;j++)
  			q_cpu[j]=-1;

  		*yy=0; *ll_cpu=ceil(d*1.0/no_of_queries_in_warp_cpu)*32;

  		cudaMalloc(&query_index_arr,sizeof(int)*d); cudaMalloc(&max_comp,sizeof(int));
  		cudaMemcpy(query_index_arr,q_cpu,sizeof(int)*d,cudaMemcpyHostToDevice);

  		cudaMalloc(&max_comp,sizeof(int));
  		cudaMemcpy(max_comp,yy,sizeof(int),cudaMemcpyHostToDevice);

      cout<<"*ll_cpu: "<<*ll_cpu<<endl;
  		cudaMalloc(&ll_gpu,sizeof(int));
  		cudaMemcpy(ll_gpu,ll_cpu,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&qq,sizeof(int));
      cudaMemcpy(qq,&q,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&d_gpu,sizeof(int));
      cudaMemcpy(d_gpu,&d,sizeof(int),cudaMemcpyHostToDevice);

      cudaMalloc(&c_gpu,sizeof(int)*d);
  		cudaMemcpy(c_gpu,c,sizeof(int)*d,cudaMemcpyHostToDevice);

      cout<<"Before srch call"<<endl;
      cout<<"ceil(d*1.0/no_of_queries_in_warp_cpu): "<<ceil(d*1.0/no_of_queries_in_warp_cpu)<<endl;
      cout<<"d: "<<d<<" ; no_of_queries_in_warp_cpu: "<<no_of_queries_in_warp_cpu<<endl;

  		srch<<<ceil(d*1.0/no_of_queries_in_warp_cpu),32>>>(c_gpu,ll_gpu,query_index_arr,max_comp,d_gpu,qq);
  		cudaDeviceSynchronize();                                    // 1 block serves no_of_queries_in_warp queries
      cout<<"After srch"<<endl;
  		cudaMemcpy(max_cmp_cpu,max_comp,sizeof(int),cudaMemcpyDeviceToHost);

  		int product=(GS_cpu)*(*max_cmp_cpu);
      sum_of_den+=product;
      printf("max_cmp_cpu: %d\n",*max_cmp_cpu);

    }
  }

}

int main()
{

  nodes=(node**)malloc(sizeof(node*)*k1);
  prefix_sum=(int*)malloc(sizeof(int)*k2);

  for(int i=0;i<k1;i++)
  {
    nodes[i]=NULL;
    if(i<k2)
      prefix_sum[i]=MAX;
  }

  int q;
  cout<<"Enter number of queries:- "<<endl;
  cin>> q;  // no_of_queries
  int** a; int *b; b=(int*)malloc(sizeof(int)*q);
  a=(int**)malloc(sizeof(int*)*q);
  for(int i=0;i<q;i++)
  {
    a[i]=(int*)malloc(sizeof(int)*(m+1));
  }
  //int a[q][m_cpu+1]; int b[q];
  for(int i=0;i<q;i++)
  {
    b[i]=MAX;
  }
  for(int i=0;i<q;i++)
  {
    for(int j=0;j<=m;j++)
      a[i][j]=MAX;
  }

  for(int i=0;i<q;i++)
  {
    cout<<"Enter 1 for inserting a tuple"<<endl;
    cout<<"Enter 2 for searching for a key"<<endl;
    cout<<"Enter 3 for doing range_query"<<endl;
    int kk; cin>>kk;
    if(kk==1)                      // insert
    {
      int dd[m]; a[i][0]=1;
      cout<<"Enter "<<m<<" elements of a tuple, of which 1st should be key:- "<<endl;
      for(int j=0;j<m;j++)
      {
        cin>>dd[j];
        a[i][j+1]=dd[j];
      }
      flag=false;
      //insert(dd);
    }
    else if(kk==2)
    {
      cout<<"Enter key to be searched"<<endl; int ke; cin>>ke;
      a[i][0]=2; a[i][1]=ke;
      //search(ke);
    }
    // else if(kk==3)
    // {
    //   cout<<"Enter 2 keys between which range_query is done"<<endl;
    //   int k1, k2; cin>>k1>>k2;
    //   range_query(k1,k2);
    // }
  }

  int jj=0;
  for(int i=0;i<q;i++)
  {
    if(a[i][0]==1)                      // b has inserts
    {
      b[jj]=i; jj++;
    }
  }

  func(a,b,q);
  if(T_count!=0)
    avg_tree_size/=T_count;

  float throughput=32*1.0/sum_of_den;
  cout<<"nk_gpu: "<<nk<<"; GS: "<<GS_cpu<<"; no_of_queries_in_warp: "<<no_of_queries_in_warp_cpu<<";  avg_tree_size: "<<avg_tree_size<<"; throughput: "<<throughput<<endl;
  //printf("throughput: %f\n",throughput);

  return 0;
}
