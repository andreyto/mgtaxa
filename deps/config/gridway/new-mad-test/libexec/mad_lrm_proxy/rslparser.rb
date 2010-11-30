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

require 'rexml/document'

# Parser of RSL2 files. It is used like this:
#
#   rsl=RSLParser.new(string_containing_the_rsl)
#   rsl.parse
#
# +rsl+ object will then have this information:
# * +executable+
# * +stdout+
# * +stderr+
# * +directory+
# * +count+
# * +jobType+
# * +environment+ (is a hash with env variable names as keys)
#
# This information can be accessed using methods with the same name:
#
#   stdout=rsl.stdout

class RSLParser
	
	# +str+ contains rsl to parse
	def initialize(str)
		@doc=REXML::Document.new(str)
		@values=Hash.new
	end
	
	# Gets data from <job><_name_></_name_></job> and adds it to
	# @values with the name converted to symbol
	def parse_value(name)
		val=@doc.elements["job/#{name.to_s}"]
		@values[name.to_sym]=val.text if val
	end
	
	# Gets all the data enclosed in <job><environment>...</environment></job>
	# and adds it as a hash to @values[:environment]
	def parse_environment
		env=Hash.new
		@doc.elements.each("job/environment") {|e|
			name=e.elements['name']
			value=e.elements['value']
			env[name.text]=value.text if name
		}
		@values[:environment]=env if env.length>0
	end
	
	# Does the parsing. It is necessary to call this method after creating the object.
	def parse
		[:executable, :stdout, :stderr, :directory, :count, 
			:jobType].each {|v| parse_value(v) }
		parse_environment
	end
	
	# Adds getters for items in @values
	def method_missing(method)
		@values[method] if @values.has_key?(method)
	end
end
