#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include <vector>
#include <random>
#include <string>
#include <climits>
#include <math.h>
#include <algorithm>
using namespace std;
unsigned test;
ifstream f ("teste.in");
ofstream g ("rezultate.txt");
string gen_message[6]= {"aleator","sortat crescator","sortat descrescator",
                        "aproape sortat crescator","aproape sortat descrescator","constant"
                       };
string sort_names[11]= {"bubble sort","count sort","merge sort","quick sort cu pivot aleator","quick sort cu pivot mediana din 3",
                        "quick sort cu pivot mediana din 5","radix sort LSD cu 128 buckets (cu matrice)","radix sort LSD cu 1024 buckets (cu matrice)",
                        "radix sort LSD cu 128 buckets","radix sort LSD cu 1024 buckets","intro sort (nativ C++)"
                       };
vector <unsigned> aux; //vector pt merge sort, din pacate un vector declarat local incetineste foarte mult
unsigned a[5];
std::mt19937 randomgen(time(NULL)); //generatorul pseudorand mersenne_twister genereaza un nr aleator in intervalul [0,10^32] si are
//mereu un seed unic.

/* tipuri de generari(gen_type):
   1:aleator
   2:sortat crescator
   3:sortat descrescator
   4:aproape sortat crescator
   5:aproape sortat descrescator
   6:constant (n_max)
*/

void generate_array(vector <unsigned> &xs,unsigned n,unsigned n_max,unsigned gen_type)
{
    switch(gen_type)
    {
        unsigned x;
    case 1:// aleator
        for (x=0; x<n; x++)
            xs[x]=randomgen()%(n_max+1); //generam n numere pseudoaleator in intervalul [0,n_max]
        break;
    case 2: //crescator
        xs[0]=randomgen()%min(n_max/10,unsigned(1000));
        for (x=1; x<n; x++)
            if (xs[x-1]<=n_max-10)
                xs[x]=xs[x-1]+randomgen()%10; //generam n numere crescatoare pseudoaleator incepand de la min(n_max/10,1000)
            else
                xs[x]=xs[x-1];
        break;
    case 3: //descrescator
        xs[0]=n_max;
        for (x=1; x<n; x++)
            if (xs[x-1]>=10)
                xs[x]=xs[x-1]-randomgen()%10; //generam n numere descrescatoare pseudoaleator inncepand cu n-max
            else
                xs[x]=xs[x-1];
        break;
    case 4: //aproape crescator
        xs[0]=randomgen()%min(n_max/10,unsigned(1000));
        for (x=1; x<n; x++)
            if (xs[x-1]<=n_max-10)
                xs[x]=xs[x-1]+randomgen()%10; //generam n numere crescatoare pseudoaleator incepand de la min(n_max/10,1000)
            else
                xs[x]=xs[x-1];
        for (x=0; x<n && n/3>0; x+=randomgen()%min(unsigned(100),n/3)+1) //il facem aproape sortat
            xs[x]=randomgen()%(n_max+1);
        break;
    case 5: //aproape descrescator
        xs[0]=n_max;
        for (x=1; x<n; x++)
            if (xs[x-1]>=10)
                xs[x]=xs[x-1]-randomgen()%10; //generam n numere descrescatoare pseudoaleator inncepand cu n-max
            else
                xs[x]=xs[x-1];
        for (x=0; x<n && n/3>0; x+=randomgen()%min(unsigned(100),n/3)+1)
            xs[x]=randomgen()%(n_max+1); //il facem aproape sortat
        break;
    case 6: //constant
        for (x=0; x<n; x++)
            xs[x]=n_max; //generam un sir constant cu valoarea n-max
    }
}

bool verif_sort(vector <unsigned> xs,unsigned n)
{
    for (unsigned x=0; x<n-1; x++)
        if (xs[x]>xs[x+1])
            return false;
    return true;
}


void bubble_sort(vector <unsigned> &xs,unsigned n)
{
    bool sorted=false;
    unsigned x;
    while(!sorted) //cat timp nu este sortat
    {
        sorted=true; //presupunem initial ca este sortat
        for (x=0; x<n-1; x++)
            if(xs[x]>xs[x+1]) //daca nu este sortat
            {
                swap(xs[x],xs[x+1]); //swap implicit c, mai rapid
                sorted=false;//sirul nu este sortat
            }
        n--; //la pasul i elementele de pe rangul [n-i+1,n] sunt pe pozitia corecta
    }
}


void count_sort(vector <unsigned> &xs,unsigned n,unsigned n_max)
{
    vector <unsigned> ap;
    unsigned x,rang=0,xmin,xmax,nr,nrap;
    xmin=xmax=xs[0];
    ap.resize(n_max+1);
    for (x=0; x<n; x++) //numaram de cate ori apare fiecare element
    {
        if (xs[x]>xmax) //retinem min si max ca sa parcurgem mai eficient
            xmax=xs[x];
        if (xs[x]<xmin)
            xmin=xs[x];
        ap[xs[x]]++;
    }
    for (nr=xmin; nr<=xmax; nr++) //construim vectorul in ordine crescatoare
        for (nrap=1; nrap<=ap[nr]; nrap++)
            xs[rang++]=nr;
}


void merge_sort (vector <unsigned> &xs,unsigned s,unsigned m,unsigned d)
{
    unsigned ps=s/*pointer stanga*/,pm=m+1/*pointer mijloc*/,rang=0,l;
    aux.resize(d-s+1);
    while (ps<=m && pm<=d) //interclasam crescator cele 2 parti ale vectorului cu ajutorul a 2 pointeri
        if(xs[ps]<xs[pm])
            aux[rang++]=xs[ps++];
        else
            aux[rang++]=xs[pm++];
    for (l=ps; l<=m; l++) //parcurgem elementele ramase
        aux[rang++]=xs[l];
    for (l=pm; l<=d; l++)
        aux[rang++]=xs[l];
    rang=0;
    for (l=s; l<=d; l++) //construim portiunea din vector ordonata
        xs[l]=aux[rang++];
}


void divide_merge(vector <unsigned> &xs,unsigned s,unsigned d)
{
    unsigned m=(s+d)/2;
    if (s<d)
    {
        divide_merge (xs,m+1,d); //impartim vectorul  in 2 parti si le sortam pe ambele cu functia mergesort
        divide_merge (xs,s,m);
        merge_sort(xs,s,m,d);
    }
}


//partitie hoare
void quick_sort(vector <unsigned> &xs,unsigned s,unsigned d,unsigned &poz_pivot)
{
    swap(xs[s],xs[poz_pivot]); //punem pivotul pe prima pozitie ca sa functioneze partitia
    unsigned pivot=xs[s];
    int ps=s-1, pd=d+1; //folosim 2 pointeri, unul la stanga altul la dreapta
    while (true)
    {
        do //gasim pozitia primului element din stanga mai mare sau egal ca pivotul
        {
            ps++;
        }
        while (xs[ps] < pivot);
        do //gasim pozitia primului element din dreapta mai mic sau egal ca pivotul
        {
            pd--;
        }
        while (xs[pd] > pivot);
        if (ps<pd) //verificam daca am aranjat dupa pivot toate elementele
            swap(xs[ps], xs[pd]); // daca nu, swap
        else
        {
            poz_pivot=pd; //altfel pivotul este acuma la pozitia pointerului din dreapta
            break;
        }
    }
}

void divide_quick(vector <unsigned> &xs,int s,int d,unsigned pivot_type)
{
    if (s<d) //while in loc de if pt partitia partilor mari
    {
        unsigned poz_pivot,m=(s+d)/2; //alegem pivotul
        switch (pivot_type)
        {
        case 1:                            //aleator
            poz_pivot=randomgen()%(d-s+1)+s;
            break;
        case 2://functiile de pseudorng mananca timp prea mult
            // mediana din 3 (aparent stanga,mijloc, dreapta, nu aleator)
            if ((xs[s] > xs[m]) xor (xs[s] > xs[d])) //aflu din 4 comparatii mediana cu xor (am gasit pe net smecheria)
                poz_pivot=s;
            else if ((xs[m] < xs[s]) xor (xs[m] < xs[d]))
                poz_pivot=m;
            else
                poz_pivot=d;
            break;
        case 3: //mediana din 5 (stanga,dreapta,mijloc,mijloc jumatatea stanga,mijloc jumatatea dreapta)
            unsigned pivot;
            //aflam mediana din max 6 comparatii fata de 10 cu sortari (alta smecherie gasita pe net)
            a[0]=s,a[1]=m,a[2]=(s+m)/2,a[3]=(m+d)/2,a[4]=d;
            if (a[0]>a[1])
                swap(a[0],a[1]);
            if (a[2]>a[3])
                swap(a[2],a[3]);
            if (a[0]<a[2])
            {
                a[0]=a[4];
                if (a[0]>a[1])
                    swap(a[0],a[1]);
            }
            else
            {
                a[2]=a[4];
                if (a[2]>a[3])
                    swap(a[2],a[3]);
            }
            if (a[0]<a[2])
            {
                if (a[1]<a[2])
                    pivot=a[2];
                else
                    pivot=a[1];
            }
            else
            {
                if (a[3]<a[0])
                    pivot=a[0];
                else
                    pivot=a[3];
            }
            if (pivot==xs[s])
                poz_pivot=s;
            else if (pivot==xs[d])
                poz_pivot=d;
            else if (pivot==xs[m])
                poz_pivot=(s+d)/2;
            else if (pivot==xs[(s+m)/2])
                poz_pivot=(s+m)/2;
            else
                poz_pivot=(m+d)/2;
            break;
        }//mediana medianelor nu mai are rost deoarece ar fi mult mai incet
        quick_sort(xs,s,d,poz_pivot);
        divide_quick(xs,s,poz_pivot,pivot_type);
        divide_quick(xs,poz_pivot+1,d,pivot_type);
    }
}


inline unsigned log2_bit(unsigned x)
{
    unsigned pow=0;
    while (x>1)
    {
        x>>=1;
        pow++;
    }
    return pow;
}

inline unsigned nr_bit_max(vector <unsigned> xs,unsigned n)
{
    unsigned vmax=xs[0],i;
    for (i=1; i<xs.size(); i++)
        if (xs[i]>vmax)
            vmax=xs[i];
    return log2_bit(vmax);

}

//e mai incet cu matrice
void radix_sort_LSD_matrix(vector <unsigned> &xs,unsigned n,unsigned base)
{
    vector <vector<unsigned>> buckets;
    unsigned slice,i,poz,bucket;
    buckets.resize(base);
    unsigned log2_base=log2_bit(base); //aflam log in baza 2 din baza (presupunem ca baza este o putere a lui 2)
    unsigned nr_max_bit=nr_bit_max(xs,n); //aflam cati biti are cel mai mare nr
    for (slice=1; slice<=nr_max_bit; slice+=log2_base) //impartim fiecare numar in sliceuri de cate log2 din baza biti,sunt introduse numerele in bucket in functie de valoare slice-ului
    {
        for (i=0; i<xs.size(); i++)
        {
            poz=((1 << log2_base) - 1) & (xs[i] >> (slice- 1));
            buckets[poz].push_back(xs[i]);
        }
        poz=0;
        for (bucket=0; bucket<buckets.size(); bucket++)
        {
            for (i=0; i<buckets[bucket].size(); i++)
                xs[poz++]=buckets[bucket][i];
            buckets[bucket].resize(0);
        }
    }
}


//nu mai folosim matrice
void radix_sort_LSD(vector <unsigned> &xs,unsigned n,unsigned base)
{
    vector <unsigned> buckets,sorted;
    unsigned slice,bucket,poz;
    int i;
    buckets.resize(base);
    sorted.resize(n);
    unsigned log2_base=log2_bit(base); //aflam log in baza 2 din baza (presupunem ca baza este o putere a lui 2)
    unsigned nr_max_bit=nr_bit_max(xs,n); //aflam cati biti are cel mai mare nr
    for (slice=1; slice<=nr_max_bit; slice+=log2_base) //impartim fiecare numar in sliceuri de cate log2 din baza biti,sunt introduse numerele in bucket in functie de valoare slice-ului
    {
        for (i=0; i<base; i++)
            buckets[i]=0;
        for (i=0; i<n; i++)
            buckets[((1 << log2_base) - 1) & (xs[i] >> (slice- 1))]++;
        for (bucket=1; bucket<buckets.size(); bucket++) //aflam ultima pozitie din vectorul sortat a ultimei cifre din galeata
            buckets[bucket]+=buckets[bucket-1];
        for (i=n-1; i>=0; i--)
        {
            poz=((1 << log2_base) - 1) & (xs[i] >> (slice- 1));
            sorted[--buckets[poz]]=xs[i]; //sortam intr-un vector auxiliar
        }
        for (i=0; i<n; i++) //mutam vectorul auxiliar sortat in vectorul original
            xs[i]=sorted[i];
    }
}

/* tipuri de sortari(sort_type)
   1:bubble sort
   2:count sort
   3:merge sort
   4:quick sort cu pivot aleator
   5:quick sort cu pivot mediana din 3
   6:quick sort cu pivot mediana din 5
   7:radix sort cu bucket 128 (cu matrice)
   8:radix sort cu bucket 1024(cu matrice)
   9:radix sort cu bucket 128
   10:radix sort cu bucket 1024
   11:intro sort
*/

void call_sort(vector <unsigned> xs,unsigned n,unsigned long long n_max,int sort_type)
{
    bool valid=true;
    g<<sort_names[sort_type-1];
    if (n_max>UINT_MAX)
    {
        sort_type=0;
        valid=false;
    }
    else if (n>100000 && sort_type==1)
        valid=false;
    else if (n_max>xs.max_size()-587000000 && sort_type==2)
        valid=false;
    else if (n>xs.max_size()-587000000)
    {
        sort_type=-1;
        valid=false;
    }
    if (valid)
    {
        auto start = std::chrono::high_resolution_clock::now();
        switch (sort_type)
        {
        case 1:
            bubble_sort(xs,n);
            break;
        case 2:
            count_sort(xs,n,n_max);
            break;
        case 3:
            divide_merge(xs,0,n-1);
            break;
        case 4:
        case 5:
        case 6:
            divide_quick(xs,0,n-1,sort_type-3);
            break;
        case 7:
            radix_sort_LSD_matrix(xs,n,128);
            break;
        case 8:
            radix_sort_LSD_matrix(xs,n,1024);
            break;
        case 9:
            radix_sort_LSD(xs,n,128);
            break;
        case 10:
            radix_sort_LSD(xs,n,1024);
            break;
        case 11:
            sort(xs.begin(),xs.end());
            break;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        g<<" sort time: "<<duration.count()<<" milisecond(s)."<<"\n";
        g<<"Memorie aditionala folosita: ";
        switch (sort_type)
        {
        case 1:
            g<<0; //Bubble sort nu foloseste memorie aditionala
            break;
        case 2:
            g<<"~"<<n_max*4; //cum folosim un vector pt a contoriza aparitiile astfel putem aproxima ca vom folosii n_max int-uri, dar count sort este O(n+k)
            break;
        case 3:
            g<<n*4; //foloseste un vector auxiliar care are maxim lungimea vectorului initial
            vector<unsigned>().swap(aux); //ne asiguram ca am eliberat zona de memorie alocata vectorului aux care e global
            //(merge doar cu swap intr-un vector temporar care o sa fie sters dupa operatie)
            break;
        case 4:
        case 5:
        case 6:
            g<<floor(log10(n)); //quick sort foloseste logn memorie aditionala
            break;
        case 7:
            g<<"~"<<4*7+n*4;
            break;
        case 8:
            g<<"~"<<4*10+n*4;
            break;
        case 9:
            g<<"~"<<4*7+n*4;
            break;
        case 10:
            g<<"~"<<4*10+n*4;
            break;
        case 11:
            g<<"? ";
            break;
        }
        g<<" byte(s)"<<"\n";
        g<<"Succes: "<<verif_sort(xs,n)<<"\n";
    }
    else
    {
        g<<" sortarea nu poate fi realizata";
        switch(sort_type)
        {
        case 0:
            g<<"(Este posibil ca numerele sa nu poata fi stocate deoarece NMax este prea mare)";
        case 1:
            g<<"(prea incet, N>10^5)";
            break;
        case 2:
            g<<"(NMax prea mare, NMax>capacitatea de stocare a clasei vector)";
            break;
        case 7:
            g<<"(N prea mare, N>capacitatea de stocare a clasei list)";
            break;
        case 8:
            g<<"(N prea mare, N>capacitatea de stocare a clasei list)";
            break;
        default:
            g<<"(N prea mare, N>capacitatea de stocare a clasei vector)";
            break;
        }
        g<<"\n";
    }
}


int main()
{
    vector <unsigned> xs;
    unsigned t,n,gen_type,sort_type;
    unsigned long long n_max;
    f>>t;
    for (test=1; test<=t; test++)
    {
        f>>n>>n_max;
        cout<<"testul "<<test<<"(N="<<n<<" NMax="<<n_max<<")"<<"\n";
        if (n<=xs.max_size()-587000000)
            xs.resize(n);
        g<<"TEST "<<test<<":"<<"\n";
        g<<"N="<<n<<"  Nmax="<<n_max<<"\n";
        for (gen_type=1; gen_type<=6; gen_type++)
        {
            cout<<"\t gen type "<<gen_message[gen_type-1]<<"\n";
            g<<"\n";
            if (n<=xs.max_size()-587000000 && n_max<=UINT_MAX)
                generate_array(xs,n,n_max,gen_type); //generam un array in functie de gen_type
            g<<"Sir "<<gen_message[gen_type-1]<<"\n\n";
            for (sort_type=1; sort_type<=11; sort_type++)
            {
                cout<<"\t \t sort type "<<sort_names[sort_type-1]<<"\n";
                call_sort(xs,n,n_max,sort_type);
                g<<"******************************"<<"\n";
            }
        }
        g<<"\n\n";
        g<<"------------------------------------";
        g<<"\n\n";
        vector<unsigned>().swap(xs); //eliberam zona de memorie a vectorului xs pt a nu avea probleme cu memoria
    }
}
