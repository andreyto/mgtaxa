# --------------------------------------------------------------------------
# Copyright 2002-2006 GridWay Team, Distributed Systems Architecture
# Group, Universidad Complutense de Madrid
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may
# not use this file except in compliance with the License. You may obtain
# a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# --------------------------------------------------------------------------

$:<<(ENV['GW_LOCATION']+'/libexec/mad_lrm_proxy')

require 'pp'
require 'gwmad'
require 'net/ssh'
require 'net/sftp'

class GWTMMadSSH < GWMad
	def initialize
		super(6,5)
		set_logger(File.open("/tmp/mad_lrm_proxy_tm.#{Process.pid}.log", "w+"))
		@threads=Array.new
	end
	
	def action_init(args)
		send_message("INIT", "-", "-", "SUCCESS")
	end
	
	def action_start(args)
		send_message("START", args[1], "-", "SUCCESS")
	end

	def action_end(args)
		send_message("END", args[1], "-", "SUCCESS")
	end
	
	def action_rmdir(args)
		(host, dir)=parse_url(args[4])
		log "#{host} #{dir}"
		
		@threads<<t=Thread.new {
			begin
				Net::SSH.start(host, ENV["USER"]) do |ssh|
					ssh.open_channel do |chan|
						chan.on_success {
							chan.close
						}
						chan.exec "rm -rf #{dir}"
					end
					ssh.loop
					ssh.close
				end
				send_message("RMDIR", args[1], "-", "SUCCESS", "-")
			rescue
				send_message("RMDIR", args[1], "-", "FAILURE", "-")
			end
			@threads.delete(t)
		}
	end

	def action_exists(args)
		send_message("EXISTS", args[1], "-", "SUCCESS")
	end
	
	def action_mkdir(args)
		(host, dir)=parse_url(args[4])
		#dir=dir[0..-2] if dir[-1].chr=='/'
		log "#{host} #{dir}"
		#begin
		@threads<<t=Thread.new {
			Net::SSH.start(host, ENV["USER"]) do |ssh|
				ssh.open_channel do |chan|
					chan.on_success {
						chan.close
					}
					chan.exec "mkdir #{dir}"
				end
				ssh.loop
				ssh.close
			end
		#rescue
			#pp e
			#send_message("RMDIR", args[1], "-", "FAILURE", "-")
		#end
			send_message("MKDIR", args[1], "-", "SUCCESS", args[4])
			@threads.delete(t)
		}
	end
	
	def action_cp(args)
		(from_host, from_dir)=parse_url(args[4])
		(to_host, to_dir)=parse_url(args[5])
		
		set_home_dir(from_dir) if from_host
		set_home_dir(to_dir) if to_host
				
		STDERR.puts "#{from_host}:#{from_dir} -> #{to_host}:#{to_dir}"
		
		session=nil
		@threads<<t=Thread.new {
		
			begin
				if from_host
					Net::SFTP.start(from_host, ENV["USER"]) do |sftp|
						session=sftp
						sftp.get_file(from_dir, to_dir)
					end
				else
					Net::SFTP.start(to_host, ENV["USER"]) do |sftp|
						session=sftp
						sftp.put_file(from_dir, to_dir)
						if args[3]=="X"
							sftp.setstat(to_dir, :permissions => 0755)
						end
					end
				end
		
				send_message("CP", args[1], args[2], "SUCCESS", "(#{args[4]}->#{args[5]})")
			rescue
				session.close if session and session.open?
				send_message("CP", args[1], args[2], "FAILURE", "(#{args[4]}->#{args[5]})")
			end
			@threads.delete(t)
		}
	end
	
	def parse_url(url)
		m=/^gsiftp\:\/\/([^\/]+)\/(.*)$/.match(url)
		if m
			[m[1], m[2]]
		elsif m=/^file\:\/\/(.*)$/.match(url)
			[nil, m[1]]
		else
			nil
		end
	end
	
	def set_home_dir(str)
		if str[0]==?~
			str[0..0]="."
		elsif str[0]!=?/
			str.insert(0, '/')
		end
		
		str
	end
	
end

mad=GWTMMadSSH.new
mad.loop
