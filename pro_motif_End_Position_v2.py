from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from sys import argv
import pandas as pd
import regex

def get_motif_seq(motif_end):#exon7's end 20bp <*****-----
    ref_file="/nfs/users/nfs_q/ql4/lustre_ql4/data/reference_sequences/CanFam3.1/Canis_familiaris.CanFam3.1.70.dna.toplevel.fa"
    ref_fa=SeqIO.parse(ref_file,"fasta")
    for seq in ref_fa:
        if seq.id=="13":
            chr13=seq.seq
            break
    motif=chr13[motif_end-1:motif_end-1+20]
    print(f"#len motif:",len(motif))
    reversed_complement_motif=motif.reverse_complement()
    print(type(motif))
    print(type(reversed_complement_motif))
    return chr13,motif,reversed_complement_motif

def get_20bp_exon(row):#different from the end or the st of intron
    exon_st=row["end"]
    exon_20bp=chr13[exon_st-20:exon_st]
    return exon_20bp

def get_other_exons_end():
    df_exon=pd.read_csv("/nfs/users/nfs_q/ql4/lustre_ql4/work/ctvt_mito/cnv_plot/new_data/annotate_genes/sample_1532T/ANGPT1_Exons.csv",sep=",",header=0)
    df_exon["nid"]=df_exon.index[::-1]+1
    df_exon["exon_start20bp"]=df_exon.apply(lambda r:get_20bp_exon(r),axis=1)
    df_exon.to_csv("ANGPT1_exon_start20bp_read.txt",sep="\t",header=True)
    return df_exon

def pro_fq_result(df,syb_name):#find split origin (only for rest read >=20bp)
    print("Nall:",len(df))
    df=df[df["nlen"]>=20]
    print("N(len of rest read>=20):",len(df))
    l_exon_syb=[]
    r1=regex.compile('(%s){s<=1,i<=1,d<=1}'%exon_X_out20bp)
    for index,row in df.iterrows():
        seq=row["rest_seq"]
        s=seq[:20]#only select the first 20bp
        print(len(s),s)
        #if s!=exon_X_out20bp:#---<*****------
        if r1.match(s)==None:
            df_exon_match=df_exon[df_exon["exon_start20bp"]==s]
            if not df_exon_match.empty:
                match_exon=df_exon_match["nid"].values[0]
                syb=match_exon
            else:
                syb="Something else!"
        else:
            syb="Unspliced"
        print(syb)
        l_exon_syb.append(syb)
    df=df.copy()
    df["splited_status"]=l_exon_syb
    print(df["splited_status"].value_counts())
    new_opt_file=f"{syb_name}.splited.ED.csv"
    df.to_csv(new_opt_file,sep=",",header=True,index=False)
    return df,new_opt_file

def extract_end1_end2_reads(fastq_file,name_list):
    df_sequences=pd.DataFrame(columns=["name","sequence"])
    l_seq_name=[]
    l_seq=[]
    #for read in SeqIO.parse(fastq_file,"fastq"):
    fastq_handle=open(fastq_file,"r")
    for title,seq,qual in FastqGeneralIterator(fastq_handle):
        seq_name=title.split("/")[0]
        if seq_name in name_list:
            l_seq_name.append(seq_name)
            l_seq.append(seq)
            print(seq)
    df_sequences["name"]=l_seq_name
    df_sequences["sequence"]=l_seq
    fastq_handle.close()
    return df_sequences

def search_bases(motif_X,seq):
    r=regex.compile('(%s){s<=2,i<=2,d<=2}'%motif_X)
    match_result=r.search(seq)
    if match_result == None:
        syb="no"# no fuzzy match"
        match_bases="N"
        sites="N"
        fuzzy_counts="N"
    else:
        syb="yes"
        match_bases=match_result.group(0)
        sites=match_result.span()
        fuzzy_counts=match_result.fuzzy_counts
    return sites,fuzzy_counts,match_bases

def search_motif_in_Seq(motif_X,seq,strand):#add strand to transfer rest bases
    #r=regex.compile('(%s){s<=2,i<=2,d<=2}'%motif_X)
    r=regex.compile('(%s){i<=1,d<=1,s<=2,2i+2d+1s<4}'%motif_X)
    match_result=r.search(seq)
    if match_result == None:
        syb="no"# no fuzzy match"
        match_bases=""
        rest_reads=""
        sites=""
        fuzzy_counts=""
    else:
        syb="yes"
        match_bases=match_result.group(0)
        sites=match_result.span()
        st=sites[0]
        ed=sites[1]
        if strand =="+":
            rest_reads=seq[:st]
        else:
            rest_reads=Seq(seq[ed:]).reverse_complement()
        fuzzy_counts=match_result.fuzzy_counts
    #print(syb)
    len_rest=len(rest_reads)
    return syb,sites,fuzzy_counts,match_bases,rest_reads,len_rest


def pro_fastq(fastq_file):
    opt_file_name=fastq_file.split("/")[-1]+".ED.txt"
    opt_file=open(opt_file_name,"w")
    opt_file.write("\t".join(["name","sequence","span","fuzzy_counts","match_bases","rest_seq","nlen","strand"])+"\n")
    nlen=len(motif)
    l_seq_end=[]
    l_seq_name=[]
    l_seq=[]
    d_seq_end_strand={}
    df_sequences=pd.DataFrame(columns=["name","sequence"])
    fastq_handle=open(fastq_file,"r")
    #for read in SeqIO.parse(fastq_file,"fastq"):
    for title,seq,qual in FastqGeneralIterator(fastq_handle):
        Seq_seq=Seq(seq)
        read_name=title
        #print(read_name)
        seq_name=read_name.split("/")[0]
        if motif in Seq_seq:
            st=Seq_seq.find(motif)
            rest_of_read=Seq_seq[:st]
            len_rest=len(rest_of_read)
            print(f"{read_name}\t+\t100%")
            opt_file.write(f"{read_name}\t{Seq_seq}\t{st}\tN\tN\t{rest_of_read}\t{len_rest}\t+\n")
            l_seq_end.append(read_name)
            d_seq_end_strand[read_name]="+"
            l_seq_name.append(seq_name)
            l_seq.append(Seq_seq)
        else:
            plus_syb,plus_sites,plus_fc,plus_mb,plus_rest_reads,len_rest1=search_motif_in_Seq(motif,seq,"+")
            if plus_syb=="yes":
                print(f"{read_name}\t+")
                opt_file.write(f"{read_name}\t{Seq_seq}\t{plus_sites}\t{plus_fc}\t{plus_mb}\t{plus_rest_reads}\t{len_rest1}\t+\n")
                l_seq_end.append(read_name)
                d_seq_end_strand[read_name]="+"
                l_seq_name.append(seq_name)
                l_seq.append(Seq_seq)
        if reversed_complement_motif in Seq_seq:
            st=Seq_seq.find(reversed_complement_motif)
            rest_of_read=Seq_seq[st+nlen:]#select all rest
            len_rest=len(rest_of_read)
            rc_rest_of_read=rest_of_read.reverse_complement()#rc it and can matchwith forward read
            print(f"{read_name}\t-\t100%")
            opt_file.write(f"{read_name}\t{Seq_seq}\t{st}\tN\tN\t{rc_rest_of_read}\t{len_rest}\t-\n")
            l_seq_end.append(read_name)
            d_seq_end_strand[read_name]="-"
            l_seq_name.append(seq_name)
            l_seq.append(Seq_seq)
        else:
            minus_syb,minus_sites,minus_fc,minus_mb,minus_rest_reads,len_rest2=search_motif_in_Seq(reversed_complement_motif,seq,"-")
            if minus_syb=="yes":
                print(f"{read_name}\t-")
                opt_file.write(f"{read_name}\t{Seq_seq}\t{minus_sites}\t{minus_fc}\t{minus_mb}\t{minus_rest_reads}\t{len_rest2}\t-\n")
                l_seq_end.append(read_name)
                d_seq_end_strand[read_name]="-"
                l_seq_name.append(seq_name)
                l_seq.append(Seq_seq)
    df_sequences["name"]=l_seq_name
    df_sequences["sequence"]=l_seq
    opt_file.close()
    fastq_handle.close()
    df_fq_result=pd.read_csv(opt_file_name,sep="\t",header=0)
    return df_sequences,l_seq_end,df_fq_result,d_seq_end_strand

def in_range_exon(st):
    df_tmp=df_exon[(df_exon["start"]<=st)&(df_exon["end"]>=st)]
    if not df_tmp.empty:
        exon=df_exon[(df_exon["start"]<=st)&(df_exon["end"]>=st)]["nid"].values[0]
    else:
        exon="Null"
    return exon

def get_paired_end_reads(dfs,names_dict):#name1_list is name_dict
    names_list=list(names_dict.keys())#list(name1_list)+list(name2_list)
    l_paired_read_st=[]
    l_split_exon=[]
    l1_mate_exon=[]
    l2_mate_exon=[]
    for index,row in dfs.iterrows():
        name1_q=row["name"]+"/1"
        name2_q=row["name"]+"/2"
        if name1_q in names_list:
            hit_strand=names_dict[name1_q]#"+"####mate_read turned out to be human read order ---->
            if hit_strand=="+":
                mate_read=Seq(row["sequence_2"]).reverse_complement()
                #sites,fuzzy_counts,match_bases=search_bases(mate_read,chr13)
            else:
                mate_read=row["sequence_2"]
            paired_read_st=chr13.find(mate_read)    
        elif name2_q in names_list:
            hit_strand=names_dict[name2_q]
            if hit_strand=="-":
                mate_read=row["sequence_1"]
            else:
                mate_read=Seq(row["sequence_1"]).reverse_complement()
            paired_read_st=chr13.find(mate_read)
        l_paired_read_st.append(paired_read_st)
        if paired_read_st!=-1:
            exon=in_range_exon(paired_read_st)
        else:
            exon="Null"
        if name1_q in names_list:
            df_mark1.loc[df_mark1["name"]==name1_q,"mate_read"]=str(mate_read)
            df_mark1.loc[df_mark1["name"]==name1_q,"mate_start"]=paired_read_st
            df_mark1.loc[df_mark1["name"]==name1_q,"mate_exon"]=exon
        elif name2_q in names_list:
            df_mark2.loc[df_mark2["name"]==name2_q,"mate_read"]=str(mate_read)
            df_mark2.loc[df_mark2["name"]==name2_q,"mate_start"]=paired_read_st
            df_mark2.loc[df_mark2["name"]==name2_q,"mate_exon"]=exon
        print(f"#paired_read_st_exon:{paired_read_st}\t{exon}")
        l_split_exon.append(exon)
    df_mark1.to_csv(file1_opt,sep=",",header=True,index=False)
    df_mark2.to_csv(file2_opt,sep=",",header=True,index=False)
    dfs["paired_start"]=l_paired_read_st
    dfs["paired_exon"]=l_split_exon
    dfs.to_csv("ANGPT_paired_read.ED.csv",sep=",",header=True,index=False)
    return dfs

if __name__=="__main__":
    exon_end=int(argv[1])#now is definitely end of exon <****-----
    fastq1_file=argv[2]
    fastq2_file=argv[3]
    chr13,motif,reversed_complement_motif=get_motif_seq(exon_end)#the end of the exon
    df_exon=get_other_exons_end()
    exon_X_out20bp=chr13[exon_end-20:exon_end]#intro adjacented to the end of exon_X
    print(f"#motif:{motif}\n#rc_motif:{reversed_complement_motif}")
    df_seq1,l_seq_end1,df_fq1_result,d_seq_end_strand1=pro_fastq(fastq1_file)
    df_seq2,l_seq_end2,df_fq2_result,d_seq_end_strand2=pro_fastq(fastq2_file)
    l_seq_end12=l_seq_end1+l_seq_end2
    l_seq_end12=[i.split("/")[0] for i in l_seq_end12]
    df_seq11=extract_end1_end2_reads(fastq1_file,l_seq_end12)    
    print(df_seq11)
    df_seq22=extract_end1_end2_reads(fastq2_file,l_seq_end12)
    print(df_seq22)
    df_sequences_end_1_2=pd.merge(df_seq11,df_seq22,on="name",suffixes=("_1","_2"))
    df_sequences_end_1_2.to_csv("detected_sequences.ED.csv",sep=",",header=True,index=False)
    df_mark1,file1_opt=pro_fq_result(df_fq1_result,"end1")
    df_mark2,file2_opt=pro_fq_result(df_fq2_result,"end2")
    df_mark1=df_mark1.copy()
    df_mark1["mate_read"]=""
    df_mark1["mate_start"]=""
    df_mark1["mate_exon"]=""

    df_mark2=df_mark2.copy()
    df_mark2["mate_read"]=""
    df_mark2["mate_start"]=""
    df_mark2["mate_exon"]=""
    names_dict={**d_seq_end_strand1,**d_seq_end_strand2}
    get_paired_end_reads(df_sequences_end_1_2,names_dict)
