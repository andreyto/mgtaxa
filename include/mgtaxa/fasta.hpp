//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the MGTAXA package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef MGT_FASTA_HPP__
#define MGT_FASTA_HPP__

namespace ADT {

	class FastaRecord {
		public:
		FastaRecord() : m_size(0) {}
		FastaRecord(const std::string& _header, 
					const std::string& _sequence = "") :
			m_head(_header), m_seq(_sequence), m_size(0) {}

		int size() const {
			int n = m_seq.size();
			if( n == 0 ) {
				n = m_size;
			}
			return n;
		}
		
		void setSize(int _size) {
			SCM_ALWAYS(m_seq.size() == 0);
			m_size = _size;
		}
		
		// The starting '>' should not be here -
		// it is eaten on read and produced on write.
		
		const std::string& getHeader() const {
			return m_head;
		}

		const std::string& getSeq() const {
			return m_seq;
		}
		
		void clear() {
			m_head = "";
			m_seq = "";
			m_size = 0;
		} 
		
		//The following members are made public only
		//for FastaReader. The client code
		//should not access them directly
		
		public:
		std::string m_head;
		std::string m_seq;
		int m_size;
		
	}; // class FastaRecord
	

	class FastaReader {
		
		public:
		
		FastaReader(const std::string& fastaFile);
		
		FastaRecord * next();
		
		int count();
		
		protected:
		
		FastaRecord m_rec;
		
		std::ifstream in;
		
		std::string m_temp;
		
		bool m_inRecord;		
		
	}; // class FastaReader

FastaReader::FastaReader(const std::string& fastaFile):
		m_inRecord(false) {
	in.open(fastaFile.c_str());
	SCM_ALWAYS(in.good());
		}

		int FastaReader::count() {
			int cnt = 0;
			for( ; this->next(); cnt++) {}
			return cnt;
		}

		FastaRecord * FastaReader::next() {

			m_rec.clear();
			m_inRecord = false;

			do {

        // If first character is >
				if (m_temp[0] == '>') {
                // If a name and a sequence were found
					if ( m_inRecord ) {
						return & m_rec;
					}
                // Sequence name isolation
					m_rec.m_head = m_temp;
					m_rec.m_head.erase(m_rec.m_head.begin());  // Character > deletion
					m_inRecord = true;
				} 
				else {
					m_rec.m_seq += m_temp;  // Sequence isolation
				}
        
			} while (std::getline(in, m_temp, '\n'));

    // Addition of the last sequence in file
			if ( m_inRecord ) {
				return & m_rec;
			}

			return 0;

		}

	
} // namespace ADT

#endif // MGT_FATSA_HPP__