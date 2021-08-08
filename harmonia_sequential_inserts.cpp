# include <bits/stdc++.h>
# define MAX 1000000000
using namespace std;

int k1=364; int k2=121;   // k1= 3^0+...+3^5 = 364= sizeof(nodes)
                          //  k2=364-3^5= sizeof(prefix_sum)
int m=2;                  // 1 key + 4 data items
int fanout=6;             // so, each node can have max 3 keys and 4 children
int nk=fanout-1;

bool flag=false;
int incremented_index=-1;

typedef struct Node
{
  int* key;      // keys are sorted
  int** data;
}node;

node** nodes;
int* prefix_sum;

int lb(int* arr, int N, int X)
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
  nn->key=(int*)malloc(sizeof(int)*nk);                // init each node with 3 keys
  nn->data=(int**)malloc(sizeof(int*)*nk);             // and 3 data arrays, each having m integers
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

int node_count=0; int prefix_sum_count=0;

int search(int key)
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

void put_in_middle(node* nn,int pos, int key, int* dd)   // put key and its dd at position pos in nn
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
  cout<<"nk: "<<nk<<endl;
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

void insert(int* dd)
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

  for(int i=0;i<q;i++)
  {
    cout<<"Enter 1 for inserting a tuple"<<endl;
    cout<<"Enter 2 for searching for a key"<<endl;
    cout<<"Enter 3 for doing range_query"<<endl;
    int kk; cin>>kk;
    if(kk==1)                      // insert
    {
      int dd[m];
      cout<<"Enter "<<m<<" elements of a tuple, of which 1st should be key:- "<<endl;
      for(int i=0;i<m;i++)
      {
        cin>>dd[i];
        //cout<<"i: "<<i<<endl;
      }
      flag=false;
      insert(dd);
    }
    else if(kk==2)
    {
      cout<<"Enter key to be searched"<<endl; int ke; cin>>ke;
      search(ke);
    }
    else if(kk==3)
    {
      cout<<"Enter 2 keys between which range_query is done"<<endl;
      int k1, k2; cin>>k1>>k2;
      range_query(k1,k2);
    }
  }

  return 0;
}
