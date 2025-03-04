##@file readsCoverage_res.py
# Compact structure to store partial results of readsCoverage for later processing; used in multi machine mode and for checkpoints.
#
#@author vlegrand@pasteur.fr
import numpy as np
import os
import time

base_chk_fname="chk_"
chk_fname_sep="_"


## Utility classes for testing the checkpoint implementation
# class checkpoint_visitor:
#     def __str__(self):
#         return self.__class__.__name__
#
# class checkpoint_visitor_11150_Cos5(checkpoint_visitor):
#     def visit(self,chk_res):
#         if chk_res.host_len!=0 or chk_res.gen!=25 or chk_res.reads_tested!=2:
#             return False
#         return True
#
# class checkpoint_visitor_38_Cos5(checkpoint_visitor):
#     def visit(self,chk_res):
#         if chk_res.host_len!=0 or chk_res.gen!=25 or chk_res.reads_tested!=2:
#             return False
#         return True






def loadArr(arr_idx0,arr_val0,arr_idx1,arr_val1,arr2D):
    for idx, val in zip(arr_idx0, arr_val0):
        arr2D[0][idx] = val

    for idx, val in zip(arr_idx1, arr_val1):
        arr2D[1][idx] = val


def loadRCRes(filename):
    npzfile = np.load(filename)
    gen_len=npzfile['gen_len']
    gen_len=int(gen_len)
    host_len=npzfile['host_len']
    host_len=int(host_len)
    termini_coverage_idx0 = npzfile['termini_coverage_idx0']
    termini_coverage_val0=npzfile['termini_coverage_val0']
    termini_coverage_idx1 = npzfile['termini_coverage_idx1']
    termini_coverage_val1 = npzfile['termini_coverage_val1']

    whole_coverage_idx0=npzfile['whole_coverage_idx0']
    whole_coverage_val0 = npzfile['whole_coverage_val0']
    whole_coverage_idx1 = npzfile['whole_coverage_idx1']
    whole_coverage_val1 = npzfile['whole_coverage_val1']

    paired_whole_coverage_idx0=npzfile['paired_whole_coverage_idx0']
    paired_whole_coverage_val0 = npzfile['paired_whole_coverage_val0']
    paired_whole_coverage_idx1 = npzfile['paired_whole_coverage_idx1']
    paired_whole_coverage_val1 = npzfile['paired_whole_coverage_val1']

    phage_hybrid_coverage_idx0=npzfile['phage_hybrid_coverage_idx0']
    phage_hybrid_coverage_val0 = npzfile['phage_hybrid_coverage_val0']
    phage_hybrid_coverage_idx1 = npzfile['phage_hybrid_coverage_idx0']
    phage_hybrid_coverage_val1 = npzfile['phage_hybrid_coverage_idx1']

    host_hybrid_coverage_idx0 = npzfile['host_hybrid_coverage_idx0']
    host_hybrid_coverage_val0 = npzfile['host_hybrid_coverage_val0']
    host_hybrid_coverage_idx1 = npzfile['host_hybrid_coverage_idx1']
    host_hybrid_coverage_val1 = npzfile['host_hybrid_coverage_val1']

    host_whole_coverage_idx0 = npzfile['host_whole_coverage_idx0']
    host_whole_coverage_val0 = npzfile['host_whole_coverage_val0']
    host_whole_coverage_idx1 = npzfile['host_whole_coverage_idx1']
    host_whole_coverage_val1 = npzfile['host_whole_coverage_val1']

    list_hybrid=npzfile['list_hybrid']
    insert=npzfile['insert']
    insert=list(insert)
    paired_mismatch=npzfile['paired_mismatch']
    reads_tested=npzfile['reads_tested']

    termini_coverage=np.array([gen_len*[0], gen_len*[0]])

    whole_coverage        = np.array([gen_len*[0], gen_len*[0]])
    paired_whole_coverage = np.array([gen_len*[0], gen_len*[0]])
    phage_hybrid_coverage = np.array([gen_len*[0], gen_len*[0]])
    host_hybrid_coverage  = np.array([host_len*[0], host_len*[0]])
    host_whole_coverage   = np.array([host_len*[0], host_len*[0]])
    loadArr(termini_coverage_idx0,termini_coverage_val0,termini_coverage_idx1,termini_coverage_val1,termini_coverage)
    loadArr(whole_coverage_idx0,whole_coverage_val0,whole_coverage_idx1,whole_coverage_val1,whole_coverage)
    loadArr(paired_whole_coverage_idx0,paired_whole_coverage_val0,paired_whole_coverage_idx1,paired_whole_coverage_val1,paired_whole_coverage)
    loadArr(phage_hybrid_coverage_idx0,phage_hybrid_coverage_val0,phage_hybrid_coverage_idx1,phage_hybrid_coverage_val1,phage_hybrid_coverage)
    loadArr(host_hybrid_coverage_idx0,host_hybrid_coverage_val0,host_hybrid_coverage_idx1,host_hybrid_coverage_val1,host_hybrid_coverage)
    loadArr(host_whole_coverage_idx0,host_whole_coverage_val0,host_whole_coverage_idx1,host_whole_coverage_val1,host_whole_coverage)

    res=RCRes(termini_coverage,whole_coverage,paired_whole_coverage,\
              phage_hybrid_coverage, host_hybrid_coverage,\
              host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested)

    return res

##
# Working structure for readsCoverage (encapsulating temporary results)
class RCWorkingS:
    def __init__(self,rc_res,cnt_line,read_match):
        self.interm_res=rc_res
        self.count_line=cnt_line
        self.read_match=read_match

class RCRes:
    def __init__(self,termini_coverage,whole_coverage,paired_whole_coverage,\
                 phage_hybrid_coverage, host_hybrid_coverage, \
                 host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested):

        self.termini_coverage=termini_coverage
        self.whole_coverage=whole_coverage
        self.paired_whole_coverage=paired_whole_coverage
        self.phage_hybrid_coverage=phage_hybrid_coverage
        self.host_hybrid_coverage=host_hybrid_coverage
        self.host_whole_coverage=host_whole_coverage

        self.list_hybrid=list_hybrid
        self.insert=insert
        self.paired_mismatch=paired_mismatch
        self.reads_tested=reads_tested

        self.gen_len = len(self.termini_coverage[0])
        self.host_len= len(self.host_hybrid_coverage[0])

    # def accept(self,a_visitor):
    #     self.vis=a_visitor
    #
    # def make_visit(self):
    #     self.vis.visit()

    def save(self,filename):
        termini_coverage_idx0 = np.flatnonzero(self.termini_coverage[0])
        termini_coverage_val0 = self.termini_coverage[0][termini_coverage_idx0]
        termini_coverage_idx1 = np.flatnonzero(self.termini_coverage[1])
        termini_coverage_val1 = self.termini_coverage[1][termini_coverage_idx1]

        whole_coverage_idx0 = np.flatnonzero(self.whole_coverage[0])
        whole_coverage_val0 = self.whole_coverage[0][whole_coverage_idx0]
        whole_coverage_idx1 = np.flatnonzero(self.whole_coverage[1])
        whole_coverage_val1 = self.whole_coverage[1][whole_coverage_idx1]

        paired_whole_coverage_idx0 = np.flatnonzero(self.paired_whole_coverage[0])
        paired_whole_coverage_val0 = self.paired_whole_coverage[0][paired_whole_coverage_idx0]
        paired_whole_coverage_idx1 = np.flatnonzero(self.paired_whole_coverage[1])
        paired_whole_coverage_val1 = self.paired_whole_coverage[1][paired_whole_coverage_idx1]

        phage_hybrid_coverage_idx0 = np.flatnonzero(self.phage_hybrid_coverage[0])
        phage_hybrid_coverage_val0 = self.phage_hybrid_coverage[0][phage_hybrid_coverage_idx0]
        phage_hybrid_coverage_idx1 = np.flatnonzero(self.phage_hybrid_coverage[1])
        phage_hybrid_coverage_val1 = self.phage_hybrid_coverage[1][phage_hybrid_coverage_idx1]

        host_hybrid_coverage_idx0 = np.flatnonzero(self.host_hybrid_coverage[0])
        host_hybrid_coverage_val0 = self.host_hybrid_coverage[0][host_hybrid_coverage_idx0]
        host_hybrid_coverage_idx1 = np.flatnonzero(self.host_hybrid_coverage[1])
        host_hybrid_coverage_val1 = self.host_hybrid_coverage[1][host_hybrid_coverage_idx1]

        host_whole_coverage_idx0 = np.flatnonzero(self.host_whole_coverage[0])
        host_whole_coverage_val0 = self.host_whole_coverage[0][host_whole_coverage_idx0]
        host_whole_coverage_idx1 = np.flatnonzero(self.host_whole_coverage[1])
        host_whole_coverage_val1 = self.host_whole_coverage[1][host_whole_coverage_idx1]

        np.savez(filename,gen_len=np.array(self.gen_len),host_len=np.array(self.host_len),\
                 termini_coverage_idx0=termini_coverage_idx0, termini_coverage_val0=termini_coverage_val0,\
                 termini_coverage_idx1=termini_coverage_idx1, termini_coverage_val1=termini_coverage_val1,\
                 whole_coverage_idx0=whole_coverage_idx0,whole_coverage_val0=whole_coverage_val0,\
                 whole_coverage_idx1=whole_coverage_idx1,whole_coverage_val1=whole_coverage_val1,\
                 paired_whole_coverage_idx0=paired_whole_coverage_idx0,paired_whole_coverage_val0=paired_whole_coverage_val0,\
                 paired_whole_coverage_idx1=paired_whole_coverage_idx1,paired_whole_coverage_val1=paired_whole_coverage_val1, \
                 phage_hybrid_coverage_idx0=phage_hybrid_coverage_idx0,phage_hybrid_coverage_val0=phage_hybrid_coverage_val0, \
                 phage_hybrid_coverage_idx1=phage_hybrid_coverage_idx1,phage_hybrid_coverage_val1=phage_hybrid_coverage_val1, \
                 host_hybrid_coverage_idx0=host_hybrid_coverage_idx0,host_hybrid_coverage_val0=host_hybrid_coverage_val0, \
                 host_hybrid_coverage_idx1=host_hybrid_coverage_idx1,host_hybrid_coverage_val1=host_hybrid_coverage_val1, \
                 host_whole_coverage_idx0=host_whole_coverage_idx0,host_whole_coverage_val0=host_whole_coverage_val0, \
                 host_whole_coverage_idx1=host_whole_coverage_idx1,host_whole_coverage_val1=host_whole_coverage_val1, \
                 list_hybrid=self.list_hybrid,insert=self.insert,paired_mismatch=np.array(self.paired_mismatch),\
                 reads_tested=self.reads_tested)


class RCCheckpoint:
    def __init__(self,count_line,core_id,idx_seq,termini_coverage,whole_coverage,paired_whole_coverage,\
                 phage_hybrid_coverage, host_hybrid_coverage, \
                 host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested,read_match):
        self.count_line=count_line
        self.core_id=core_id
        self.idx_seq=idx_seq
        self.read_match=read_match
        self.res=RCRes(termini_coverage,whole_coverage,paired_whole_coverage,\
                 phage_hybrid_coverage, host_hybrid_coverage, \
                 host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested)


    def save(self,dir_chk,core_id,idx_refseq):
        filename=base_chk_fname+str(self.core_id)+chk_fname_sep+str(self.idx_seq)+chk_fname_sep+\
                 str(self.count_line)+chk_fname_sep+str(self.read_match)
        full_fname = os.path.join(dir_chk, filename)
        self.res.save(full_fname)
        # remove old breakpoint file
        list_f=os.listdir(dir_chk)
        sub_s=base_chk_fname+ str(core_id) + chk_fname_sep + str(idx_refseq) + chk_fname_sep
        for f in list_f:
            if f!=filename+".npz" and sub_s in f:
                os.remove(os.path.join(dir_chk,f))


class RCCheckpoint_handler:
    def __init__(self,chk_freq,dir_chk,test_mode=False):
        self.chk_freq=chk_freq
        self.test_mode = test_mode
        self.start_t=0
        self.dir_chk = dir_chk
        # if self.test_mode == True:
        #     self.v38_C5 = checkpoint_visitor_38_Cos5()
        #     self.v11150_C5 = checkpoint_visitor_11150_Cos5()
        if self.test_mode==True:
            self.start_t = time.perf_counter_ns()
            if os.path.exists(dir_chk):
                if not (os.path.isdir(dir_chk)):
                    raise RuntimeError("dir_chk must point to a directory")
            else:
                os.mkdir(dir_chk)
        elif self.chk_freq!=0:
            if os.path.exists(dir_chk):
                if not (os.path.isdir(dir_chk)):
                    raise RuntimeError("dir_chk must point to a directory")
            else:
                raise RuntimeError("dir_chk must point to an existing directory")

    def getIdxSeq(self,core_id):
        idx_seq=0
        if self.chk_freq!=0 or self.test_mode==True:
            list_f = os.listdir(self.dir_chk)
            subfname = base_chk_fname+ str(core_id) + chk_fname_sep
            chk_f = ""
            for fname in list_f:
                if subfname in fname:
                    chk_f = fname
                    break
            if chk_f != "":
                l=chk_f.split(chk_fname_sep)
                idx_seq=int(l[2])
        return idx_seq


    def load(self,core_id,idx_refseq):
        if self.chk_freq!=0 or self.test_mode==True:
            list_f = os.listdir(self.dir_chk)
            subfname = base_chk_fname+ str(core_id) + chk_fname_sep + str(idx_refseq) + chk_fname_sep
            chk_f = ""
            for fname in list_f:
                if subfname in fname:
                    chk_f = fname
                    break
            if chk_f != "":
                interm_res=loadRCRes(os.path.join(self.dir_chk,chk_f))
                # if self.test_mode==True:
                #     interm_res.accept(self.v38_C5)
                l=chk_f.split(chk_fname_sep)
                cnt_line=int(l[-2])
                tmp=l[-1] # get rid of .npz extension
                l2=tmp.split(".")
                read_match=int(l2[0])
                partial_res=RCWorkingS(interm_res,cnt_line,read_match)
                # if self.test_mode:
                #     partial_res.accept(self.v38_C5)
                #     partial_res.make_visit()
                return partial_res
            else:  # no checkpoint found for this sequence, start from beginning
                return None
        else:
            return None


    def check(self,count_line,core_id,idx_seq,termini_coverage,whole_coverage,paired_whole_coverage,\
                 phage_hybrid_coverage, host_hybrid_coverage, \
                 host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested,read_match):
        cur_t = time.perf_counter_ns()
        elapsed_t = cur_t - self.start_t
        #convert elapsed_t tp to seconds
        elaspsed_t=elapsed_t * 0.000000001
        if (self.test_mode==True or (self.chk_freq!=0 and (elapsed_t % self.chk_freq == 0))):  # time to create checkpoint.
            chkp=RCCheckpoint(count_line,core_id,idx_seq,termini_coverage,whole_coverage,paired_whole_coverage,\
                 phage_hybrid_coverage, host_hybrid_coverage, \
                 host_whole_coverage,list_hybrid,insert,paired_mismatch,reads_tested,read_match)
            chkp.save(self.dir_chk,core_id,idx_seq)


    def end(self,core_id):
        if (self.test_mode==False and self.chk_freq!=0) :
            # remove old breakpoint file
            list_f = os.listdir(self.dir_chk)
            sub_s=base_chk_fname+str(core_id)+chk_fname_sep
            for f in list_f:
                if sub_s in f:
                    os.remove(os.path.join(self.dir_chk, f))









