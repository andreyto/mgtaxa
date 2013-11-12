#!/usr/bin/env python
from subprocess import check_call
import sh
import sys, os, glob, shutil, tempfile
from contextlib import nested
from os.path import join as pjoin
import ftputil, fnmatch
from MGT.Common import *
from MGT.BatchMakeflow import *
from MGT.Logging import *

from argh import ArghParser,arg

log = logging.getLogger(os.path.basename(sys.argv[0]))

class ftputil_download_progress_reporter:

    import time

    def __init__(self,host,remote_file,local_file,log=None,log_period=10,keep_alive=False):
        self.host = host
        self.remote_file = remote_file
        self.local_file = local_file
        self.log = log
        self.log_period = log_period
        self.keep_alive = keep_alive
        self.length = 0
        self._last_cnt = 0
        self._last_sec = time.time()

    def __call__(self,chunk):
        self.length += len(chunk)
        if self.log:
            cnt = int(self.length/2**20)
            if cnt > self._last_cnt:
                self._last_cnt = cnt
                sec = time.time()
                if sec > self._last_sec + self.log_period:
                    self._last_sec = sec
                    self.log.info("{}: {}MB ready".\
                        format(self.local_file,cnt))
        if keep_alive:
            self.host.keep_alive()

class RefseqDownloader:

    def __init__(self):
        self.wrapper = options["wrapper"]

    def list_remote_files(self,
            config_file,
            remote_paths_file
            ):
        conf = load_config_json(config_file)
        conf_ftp = conf["ftp"]
        conf_db = conf_ftp["refseq"]
        path_base = conf_db["path_base"]
        file_paths = []
        with ftputil.FTPHost(conf_ftp["host"], 'anonymous', options.toolEmail) as host:
            for org_group in conf_db["org_groups"]:
                path_group = host.path.join(path_base,org_group)
                host.chdir(path_group)
                names = host.listdir(host.curdir)
                for file_type in conf_db["file_types"]:
                    patt = "*.{}.gz".format(file_type)
                    names_type = fnmatch.filter(names,patt)
                    for name_type in names_type:
                        do_excl = False
                        for exclude_rx in conf_db["exclude_rxs"]:
                            if re.match(exclude_rx,name_type):
                                do_excl = True
                                break
                        if not do_excl:
                            file_paths.append(dict(
                                org_group=org_group,
                                file_type=file_type,
                                r_file=host.path.join(path_group,name_type)
                                ))
        save_config_json(file_paths,remote_paths_file)

    def wf_download_files(self,
            config_file,
            remote_paths_file,
            local_paths_file,
            makeflow_file,
            output_dir
            ):
        conf = load_config_json(config_file)
        conf_ftp = conf["ftp"]
        wrapper = self.wrapper
        file_paths = load_config_json(remote_paths_file)
        path_hasher = PathHasher(output_dir)
        mkf = MakeflowWriter(makeflow_file,mode="a")
        cmd_pref = cmd_self(wrapper)["cmd"]
        for rec in file_paths:
            local_file = path_hasher(rec["r_file"].split("/")[-1])
            cmd = "{} internal download-file {} {} {}".\
                    format(cmd_pref,config_file,
                    rec["r_file"],local_file)
            mkf.appendJob(
                    targets=[local_file],
                    inputs=[config_file],
                    cmd=cmd,
                    vars=[
                        "@BATCH_LOCAL=1"
                        ]
                    )
            rec["l_file"] = local_file
        save_config_json(file_paths,local_paths_file)
    
    def download_file(self,
            config_file,
            remote_file,
            local_file
            ):
        conf = load_config_json(config_file)
        conf_ftp = conf["ftp"]

        with ftputil.FTPHost(conf_ftp["host"], 'anonymous', options.toolEmail) as host:
            host.download(remote_file,local_file,mode="b",
                    callback=ftputil_download_progress_reporter(
                        host=host,
                        remote_file=remote_file,
                        local_file=local_file,
                        log=log)
                    )
    
    def extract_fasta_header(self,fasta_file,hdr_file):
        inp = openCompressed(fasta_file,"r")
        out = openCompressed(hdr_file,"w")
        for line in inp:
            if line.startswith(">"):
                out.write(line)
        inp.close()
        out.close()

    def wf_extract_headers(self,
            config_file,
            local_paths_file,
            hdr_paths_file,
            makeflow_file
            ):
        conf = load_config_json(config_file)
        wrapper = self.wrapper
        file_paths = load_config_json(local_paths_file)
        mkf = MakeflowWriter(makeflow_file,mode="a")
        cmd_pref = cmd_python(wrapper)["cmd"]
        cmd_self_pref = cmd_self(wrapper)["cmd"]
        for rec in file_paths:
            local_file = rec["l_file"]
            do_skip = False
            if fnmatch.fnmatch(local_file,"*.g[bp]ff.gz"):
                hdr_file = stripSfx(local_file)+".hdr.gz"
                cmd = "{} $MGT_HOME/bin/extractGbHeader.py {} {}".\
                        format(cmd_pref,local_file,hdr_file)
            elif fnmatch.fnmatch(local_file,"*.f[na]a.gz"):
                hdr_file = stripSfx(local_file)+".hdr.gz"
                cmd = "{} internal extract-fasta-header {} {}".\
                        format(cmd_self_pref,local_file,hdr_file)
            else:
                do_skip = True
            if not do_skip:
                mkf.appendJob(
                        targets=[hdr_file],
                        inputs=[config_file,local_file],
                        cmd=cmd
                        )
                rec["hdr_file"] = local_file
            else:
                rec["hdr_file"] = None
        save_config_json(file_paths,hdr_paths_file)

    def wf_generate(self,
            config_file,
            remote_paths_file,
            local_paths_file,
            hdr_paths_file,
            makeflow_file,
            output_dir
            ):
        """Generate internal workflow"""
        #create default header
        MakeflowWriter(makeflow_file).close()
        #these will append to the makeflow file
        self.wf_download_files(
            config_file,
            remote_paths_file,
            local_paths_file,
            makeflow_file,
            output_dir
            )
        self.wf_extract_headers(
            config_file,
            local_paths_file,
            hdr_paths_file,
            makeflow_file
            )

    @arg("config-file",type=os.path.abspath)
    @arg("output-dir",type=os.path.abspath)
    def checkout(self,
            config_file,
            output_dir
            ):
        """Main entry point: get and process the database"""
        conf = load_config_json(config_file)
        wrapper = self.wrapper
        cmd = cmd_self(wrapper)
        cmd_pref = cmd["cmd"]
        name_root = stripSfx(os.path.basename(cmd["script"]))

        remote_paths_file = name_root+".remote_paths.json"
        local_paths_file = name_root+".local_paths.json"
        hdr_paths_file = name_root+".hdr_paths.json"
        makeflow_file = name_root+".internal.mkf"
        makeflow_file_checkout = name_root+".mkf"

        mkf = MakeflowWriter(makeflow_file_checkout)
        mkf.appendJob(
                targets=[remote_paths_file],
                inputs=[config_file],
                cmd="{} internal list-remote-files {} {}".format( 
                    cmd_pref,
                    config_file,
                    remote_paths_file
                    )
                )
        mkf.appendJob(
                targets=[
                    local_paths_file,
                    hdr_paths_file,
                    makeflow_file
                    ],
                inputs=[
                    config_file,
                    remote_paths_file
                    ],
                cmd="{} internal wf-generate {} {} {} {} {}".format( 
                    cmd_pref,
                    config_file,
                    remote_paths_file,
                    local_paths_file,
                    hdr_paths_file,
                    makeflow_file,
                    output_dir
                    )
                )

        mkf.appendMakeflow(
                    flow = makeflow_file
                    )
        mkf.close()


def main():
    logging_config(detail="high")
    parser = ArghParser()
    ref_down = RefseqDownloader()
    parser.add_commands([
        ref_down.checkout
        ]
        )
    parser.add_commands([
        ref_down.list_remote_files,
        ref_down.download_file,
        ref_down.wf_generate,
        ref_down.extract_fasta_header
        ],
        namespace="internal",
        title="Commands used internally by the application"
        )
    parser.dispatch()

if __name__ == "__main__":
    main()

