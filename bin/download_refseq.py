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

log = logging.getLogger(os.path.basename(sys.argv[0]))

class ftputil_download_progress_reporter:

    import time

    def __init__(self,host,remote_file,local_file,log=None,log_period=10):
        self.host = host
        self.remote_file = remote_file
        self.local_file = local_file
        self.log = log
        self.log_period = log_period
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
        self.host.keep_alive()

class RefseqDownloader:

    def list_files(self,
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
            makeflow_file
            ):
        conf = load_config_json(config_file)
        conf_ftp = conf["ftp"]
        wrapper = conf.get("wrapper","")
        file_paths = load_config_json(remote_paths_file)
        mkf = MakeflowWriter(makeflow_file,mode="a")
        cmd_pref = cmd_self(wrapper)["cmd"]
        for rec in file_paths:
            local_file = rec["r_file"].split("/")[-1]
            cmd = "{} download-file {} {} {}".\
                    format(cmd_pref,config_file,
                    rec["r_file"],local_file)
            mkf.appendJob(
                    targets=[local_file],
                    inputs=[config_file],
                    cmd=cmd
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

    def wf_extract_headers(self,
            config_file,
            local_paths_file,
            hdr_paths_file,
            makeflow_file
            ):
        conf = load_config_json(config_file)
        wrapper = conf.get("wrapper","")
        file_paths = load_config_json(local_paths_file)
        mkf = MakeflowWriter(makeflow_file,mode="a")
        cmd_pref = cmd_python(wrapper)["cmd"]
        for rec in file_paths:
            local_file = rec["l_file"]
            do_skip = False
            if fnmatch.fnmatch(local_file,"*.g[bp]ff.gz"):
                hdr_file = stripSfx(local_file)+".hdr.gz"
                cmd = "{} $MGT_HOME/bin/extractGbHeader.py {} {}".\
                        format(cmd_pref,local_file,hdr_file)
            elif fnmatch.fnmatch(local_file,"*.f[na]a.gz"):
                hdr_file = stripSfx(local_file)+".hdr.gz"
                cmd = "{} zgrep '>' {} | gzip -c > {}".\
                        format(wrapper,local_file,hdr_file)
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

def main():
    from argh import ArghParser
    logging_config(detail="high")
    parser = ArghParser()
    ref_down = RefseqDownloader()
    parser.add_commands([
        ref_down.list_files,
        ref_down.download_file,
        ref_down.wf_download_files,
        ref_down.wf_extract_headers
        ])
    parser.dispatch()

if __name__ == "__main__":
    main()

