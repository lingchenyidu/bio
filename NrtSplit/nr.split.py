from Bio import SeqIO
import pandas as pd
import sys
class nr_split(object):
    def __init__(self):
        self.nodes_names = ['taxid','num1','num2','num3','num4','num5','num6','num7','num8','num9','num10',
               'num11','num12','num13']
    def read_nodes(self,nodesfile):
        nodes_data = pd.read_csv(nodesfile,sep='\t\|\t',names=self.nodes_names,engine='python')
        nodes_data['id'] = nodes_data['num4']
        return nodes_data
    def read_division(self,divisionfile,nodes_data):
        division_data = pd.read_csv(divisionfile,sep='\t\|\t',names=['id','sx','name','other'])
        division_data_tmp = nodes_data.merge(division_data,on='id')[['taxid','sx','name']]
        return division_data_tmp
    def get_asskey(self,protfile,divisionfile,nodesfile):
        nodes_data = self.read_nodes(nodesfile)
        division_data_tmp = self.read_division(divisionfile,nodes_data)
        alldata = pd.DataFrame({})
        data = pd.read_table(protfile,iterator=True)
        while True:
            try:
                tmpdata = data.get_chunk(1000000)
                alltmpdata = tmpdata.merge(division_data_tmp,on='taxid')[['accession.version','sx']]
                alldata = pd.concat([alldata,alltmpdata])
                # break
            except StopIteration:
                print("Iteration is stopped.")
                break
        return alldata.set_index('accession.version').to_dict()['sx']
    def get_fa(self,nrfile,asskeydict,outdir):
        for seq_record in SeqIO.parse(nrfile, "fasta"):
            accid = str(seq_record.id)
            if accid in asskeydict:
                filesave = open(outdir+'/'+asskeydict[accid]+'.fa','a+')
                filesave.write('>'+str(seq_record.description)+'\n')
                filesave.write(str(seq_record.seq)+'\n')
                filesave.close()
            else:
                filesave = open(outdir+'/unkonwn.fa','a+')
                filesave.write('>'+str(seq_record.description)+'\n')
                filesave.write(str(seq_record.seq)+'\n')
                filesave.close()
    def main(self):
        import argparse,sys
        print(sys.argv)
        if len(sys.argv) < 2:
            print('please python3 ./nr.split.py --help')
            sys.exit(1)
        usage = 'nr.split.py [--outdir] [--nrfile] [--nodesfile] [--divisionfile] [--protfile]'
        p = argparse.ArgumentParser(usage='python3 ./nr.split.py [--outdir] [--nrfile] [--nodesfile] [--divisionfile] [--protfile]', description='')
        p.add_argument('-o', '--outdir', type=str, help='outdir')
        p.add_argument('-n', '--nrfile', type=str, help='nrfile or ntfile')
        p.add_argument('-t', '--nodesfile', type=str, help='nodes file')
        p.add_argument('-d', '--divisionfile', type=str, help='division  file')
        p.add_argument('-p', '--protfile', type=str, help='prot file')
        args = p.parse_args()
        outdir = args.outdir
        nrfile = args.nrfile
        nodesfile = args.nodesfile
        divisionfile = args.divisionfile
        protfile = args.protfile
        asskeydict = self.get_asskey(protfile,divisionfile,nodesfile)
        self.get_fa(nrfile,asskeydict,outdir)

if __name__ == "__main__":
    Run = nr_split()
    run = Run.main()