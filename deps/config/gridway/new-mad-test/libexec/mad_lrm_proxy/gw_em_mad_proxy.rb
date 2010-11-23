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

require 'gwmad'
require 'rslparser'
require 'sshresource'
require 'sshjob'

class GWEMMadSSH < GWMad
	def initialize
		super(4)
		set_logger(File.open("/tmp/mad_lrm_proxy.#{Process.pid}.log", "w"))
		@jobs=Hash.new
		@resources=Hash.new
	end
	
	def action_init(args)
		send_message("INIT", "-", "SUCCESS")
	end
	
	def action_submit(args)
		parsed_job=RSLParser.new(File.new(args[3]))
		parsed_job.parse
		
		host=args[2].split('/')[0]
		user=ENV['USER']
		
		job=SSHJob.new(
					parsed_job.executable,
					"",
					parsed_job.stdout,
					parsed_job.stderr
		)
		
		job.environment=parsed_job.environment
		job.gw_id=args[1]
		
		# Save the job to query its state
		@jobs[args[1]]=job
		
		resource=nil
		if @resources.has_key?(host)
			resource=@resources[host]
		else
			resource=SSHResource.new(host, user)
			@resources[host]=resource
		end
		
		resource.start_loop
		resource.submit(job)
		
		send_message("SUBMIT", args[1], "SUCCESS", "ssh://#{user}@#{host}/#{job.hash}")
	end
	
	def action_poll(args)
		if @jobs.has_key?(args[1])
			send_message("POLL", args[1], "SUCCESS", @jobs[args[1]].status)
		end
	end
end

mad=GWEMMadSSH.new
mad.loop
