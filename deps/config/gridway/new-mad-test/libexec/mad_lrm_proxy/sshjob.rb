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

require 'digest/md5'

class SSHJob
	attr_reader :hash
	attr_accessor :status, :time_start, :time_end, :exit_code, :gw_id
	attr_accessor :environment
	
	def initialize(executable, parameters, stdout, stderr)
		@executable=executable
		@parameters=parameters
		@stdout=stdout
		@stderr=stderr
		@environment=nil
		@status='PENDING'
		@hash=gen_hash
		@remote_wrapper="/tmp/ssh_wrapper.#{@hash}"
		
		@script_vars={
			'__IDENTIFIER__' => @hash,
			'__STDOUT__' => @stdout,
			'__STDERR__' => @stderr,
			'__COMMAND__' => "#{@executable} #{@parameters}",
			'__ENVIRONMENT__' => ""
		}
	end
	
	def submit(session)
		@session=session
					
		line="nohup "
		line+="#{@remote_wrapper} #{@hash} #{@stdout} #{@stderr} "
		line+=@executable+" "
		line+=@parameters+" " if @parameters
		line+="2>/tmp/ssh_wrap.#{@hash}.err >/tmp/ssh_wrap.#{@hash}.out &"
		
		@session.open_channel {|chan|
		
			chan.on_success {
				chan.send_data "cat <<EOT >#{@remote_wrapper}\n"
				
				script=""
				File.open(ENV["GW_LOCATION"]+"/libexec/ssh_wrapper.sh", "r") do |f|
					script=f.read
				end
				
				if @environment
					env=""
					@environment.each {|key, val|
						env+="export #{key}=\"#{val}\"\n"
					}
					@script_vars['__ENVIRONMENT__']=env
				end
				
				@script_vars.each {|key, val|
					script.gsub!(key, val)
				}
				
				chan.send_data script
				
				chan.send_data "\nEOT\n"
				chan.send_data "chmod +x #{@remote_wrapper}\n"
			
				chan.send_data "#{line}\n"
				chan.send_data "exit\n"
				chan.close
				@session.close
			}
		
			chan.send_request "shell", nil, true		
		}
		
		while @session.open?
			session.connection.process false
		end	
	end
	
	def send_wrapper
		puts "sending #{@remote_wrapper}"
		@session.sftp.connect do |sftp|
			sftp.put_file 'ssh_wrapper.sh', @remote_wrapper
			sftp.setstat(@remote_wrapper, :permissions => 0700)
		end
		puts "sent #{@remote_wrapper}"
	end
	
	def delete_wrapper
		@session.sftp.connect do |sftp|
			sftp.remove(@remote_wrapper)
		end
	end
	
	def gen_hash
		# Generate actual time in epoch format and convert it to string
		time=Time.now.to_i.to_s
		Digest::MD5.hexdigest(
			"#{time} #{@executable} #{@parameters} #{@stdout} #{@stderr} #{rand}")
	end
		
end
