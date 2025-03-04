##@file debug_utils.py
#
# Contains utility classes for debugging and testing.
#
#@author vegrand@pasteur.fr

## Utility class for debugging.
#
# Contains the mapping results for 1 read.
# For map_start,map_end,map_rcpl_start,map_rcpl_stop a value of 0 means that no match was found and a value of 1 means that a match was found.
class ReadMappingInfo:
    ##
    #
    # @param idx_read Number of the read in the processing (reads are processed in the same order as they are found in the fasta file).
    # @param map_start Read maps at its beginning (1rts seed characters) or not.
    # @param map_end Read maps at is end (last seed characters) or not.
    # @param map_rcpl_start Start of reverse complement maps or not.
    # @param map_rcpl_end End of reverse complement maps or not.
    def __init__(self,idx_read,map_start,map_end,map_rcpl_start,map_rcpl_stop):
        self.idx_read=idx_read
        self.map_start=map_start
        self.map_end=map_end
        self.map_rcpl_start=map_rcpl_start
        self.map_rcpl_end=map_rcpl_stop

        pass




## Aim of this class is to give the ability to compare the results of readsCoverage (oriinal CPU version) and readsCoverageGPU.
class ReadMappingInfoLogger:
    ##
    #
    # @param cnt_read count only reads that were not rejected (readlen >= seed)
    # @param l_rm_info list of ReadMappingInfo objects.
    # @param cur_r_info ReadMappingInfo for read that s currently being processed.
    def __init__(self):
        self.cnt_read = 0 # count only reads that were not rejected (readlen >= seed)
        self.l_rm_info=[]
        self.cur_r_info=None
        self.rw_lst = []

    def add_rw(self, rw):
        self.rw_lst.append(rw)

    def newRmInfo(self,numR_in_file=None):
        if self.cur_r_info!=None:
            self.l_rm_info.append(self.cur_r_info)
        if (numR_in_file!=None):
            idx_read=numR_in_file
        else:
            idx_read=self.cnt_read
        self.cur_r_info=ReadMappingInfo(idx_read,0,0,0,0)
        self.cnt_read+=1

    ## Records the mapping information (does it map or not and where) for the read that is currently being processed.
    def rMatch(self,akey):
        if self.cur_r_info == None:
            raise RuntimeError("Call newRmInfo() before calling rsMatch()")
        if akey=="mstart":
            self.cur_r_info.map_start = 1
        elif akey=="mend":
            self.cur_r_info.map_end=1
        elif akey=="mrcplstart":
            self.cur_r_info.map_rcpl_start=1
        elif akey=="mrcplend":
            self.cur_r_info.map_rcpl_end=1
        else:
            raise RuntimeError("invalid key to indicate where read matches sequence")

    def getMatchInfoList(self):
        if self.cur_r_info != None:
            self.l_rm_info.append(self.cur_r_info)
        return self.l_rm_info

    ## Flushes all ReadMappingInfo to the given file.
    def flush(self,filename):
        self.f_debug = open(filename, "w")
        if self.cur_r_info != None:
            self.l_rm_info.append(self.cur_r_info)
        self.f_debug.write(self.cnt_read)
        for elm in self.l_rm_info:
            my_str=str(elm.idx_read)+"|"+str(elm.map_start)+"|"+str(elm.map_end)+"|"+str(elm.map_rcpl_start)+"|"+str(elm.map_rcpl_stop)
            self.f_debug.write(my_str)
        self.f_debug.close()

