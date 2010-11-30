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

require 'net/ssh'
require 'net/sftp'

class SSHResource
	def initialize(hostname, username)
		@hostname=hostname
		@username=username
		@session=nil
		@jobs=Hash.new
		@threads=Array.new
		@loop_thread=nil
		start_job_monitor		
	end
	
	def submit(job)
		@jobs[job.hash]=job
		@threads<<Thread.new {
			job.submit(open_session)
		}
	end
	
	def session
		@session=open_session if !@session or !@session.open?
		@session
	end
	
	def open_session
		Net::SSH.start(@hostname, @username)
	end

	def process_monitor_line(line)
		(time, hash, state, exit_code)=line.strip.split(':')
		STDERR.puts line
		#STDERR.puts @jobs.to_s
		if @jobs.has_key?(hash)
			job=@jobs[hash]
			if state=='start'
				job.time_start=time
				job.status='ACTIVE'
			elsif state=='end'
				job.time_end=time
				job.exit_code=exit_code
				if exit_code=="0"
					job.status='DONE'
				else
					job.status='FAILED'
				end
			end
			STDERR.puts "#{hash} -> #{job.status}"
			
			# Print CALLBACK
			# TODO: Move this to mad code and add custom callbacks
			STDOUT.puts "CALLBACK #{job.gw_id} SUCCESS #{job.status}"
			STDERR.puts "CALLBACK #{job.gw_id} SUCCESS #{job.status}"
			STDOUT.flush
			STDERR.flush
		end
		STDERR.flush
	end
	
	def start_job_monitor
		session.open_channel {|chan|
			chan.on_data {|c, data|
				data.each_line {|line|
					#STDERR.puts "  Monitor Data: #{line}"
					self.process_monitor_line(line)				
				}
			}
			chan.exec("touch ~/.gw_ssh_status && " <<
				"tail -f ~/.gw_ssh_status")
		}
	end
	
	def start_loop
		@loop_thread=Thread.new do
			while true
				session.connection.process false
				STDERR.puts "LOOP!"
			end
		end if !@loop_thread
	end

end
