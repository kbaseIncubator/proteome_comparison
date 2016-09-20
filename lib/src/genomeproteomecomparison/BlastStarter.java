package genomeproteomecomparison;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

public class BlastStarter {
	private static final Pattern tabDiv = Pattern.compile(Pattern.quote("\t"));
	private static long lastTimestamp = -1;
	
	public static void run(File tempDir, Map<String, String> proteomeQ, Map<String, String> proteome2, 
			File blastBin, String maxEvalue, ResultCallback ret) throws Exception {
		File q_file = saveFastaSeries(tempDir, proteomeQ, "query");
		File s_file = null;
		if (proteome2 != null) {
			s_file = saveFastaSeries(tempDir, proteome2, "db");
		} else {
			s_file = q_file;
		}
		try {
			runBlast(blastBin, q_file, s_file, tempDir, maxEvalue, ret);
		} finally {
			Thread.sleep(200);
			q_file.delete();
			if (proteome2 != null) {
				s_file.delete();
			}
			File s2_file = new File(tempDir, s_file.getName()+".phr");
			if (s2_file.exists()) 
				s2_file.delete();
			File s3_file = new File(tempDir, s_file.getName()+".pin");
			if (s3_file.exists()) 
				s3_file.delete();
			File s4_file = new File(tempDir, s_file.getName()+".psq");
			if (s4_file.exists()) 
				s4_file.delete();
		}
	}

	private synchronized static long generateTimestmap() {
		long ret = System.currentTimeMillis();
		if (ret <= lastTimestamp)
			ret = lastTimestamp + 1;
		lastTimestamp = ret;
		return ret;
	}
	
	private static File saveFastaSeries(File tempDir, Map<String, String> proteome, String prefix) throws Exception {
		File ret = new File(tempDir, "tmp_"+prefix+"_"+generateTimestmap()+".fst");
		FastaWriter fw = new FastaWriter(new PrintWriter(ret));
		for (Map.Entry<String, String> entry : proteome.entrySet())
			fw.write(entry.getKey(), entry.getValue());
		fw.close();
		return ret;
	}
	
	private static void runBlast(File blastBin, File queryFile, File databaseFile, 
			File tempDir, String maxEvalue, ResultCallback ret) throws Exception {
		File tempResFile = new File(tempDir, "tmp_result_" + generateTimestmap() + ".txt");
		CorrectProcess cp = null;
		ByteArrayOutputStream err_baos = null;
		Exception err = null;
		String binPath = blastBin == null ? "" : (blastBin.getAbsolutePath() + "/");
		try {
			Process p = Runtime.getRuntime().exec(CorrectProcess.arr(binPath + "makeblastdb", 
					"-dbtype", "prot", "-in", databaseFile.getAbsolutePath()));
			err_baos = new ByteArrayOutputStream();
			cp = new CorrectProcess(p,new ByteArrayOutputStream(),"",err_baos,"");
			p.waitFor();
			err_baos.close();
		}catch(Exception ex) {
			try{ err_baos.close(); }catch(Exception ex_) {}
			try{ if(cp!=null)cp.destroy(); }catch(Exception ex_) {}
			err = ex;
		}
		if(err_baos!=null) {
			String err_text = new String(err_baos.toByteArray());
			if(err_text.length()>0)
				err = new Exception("makeblastdb: "+err_text,err);
		}
		int procExitValue = -1;
		if(err==null) {
			cp = null;
			err_baos = null;
			try {
				Process p = Runtime.getRuntime().exec(CorrectProcess.arr(binPath + "blastp",
						"-query", queryFile.getAbsolutePath(), "-db", databaseFile.getAbsolutePath(), 
						"-out", tempResFile.getAbsolutePath(), "-outfmt", "6", "-evalue", maxEvalue));
				err_baos = new ByteArrayOutputStream();
				cp = new CorrectProcess(p,new ByteArrayOutputStream(),"",err_baos,"");
				p.waitFor();
				err_baos.close();
				procExitValue = p.exitValue();
			}catch(Exception ex) {
				try{ err_baos.close(); }catch(Exception ex_) {}
				try{ if(cp!=null)cp.destroy(); }catch(Exception ex_) {}
				err = ex;
			}
			if(err_baos!=null) {
				String err_text = new String(err_baos.toByteArray());
				if(err_text.length()>0)
					err = new Exception("blastp: "+err_text,err);
			}
		}
		if (procExitValue == 0 && tempResFile.exists()) {
			loadData(tempResFile,ret);
			tempResFile.delete();		
		} else {
			if (err == null)
				err = new IllegalStateException("Blast exit code: " + procExitValue + ", " +
						"result file existance: " + tempResFile.exists());
			if (tempResFile.exists()) 
				tempResFile.delete();		
			throw err;
		}
	}
	
	private static void loadData(File f, ResultCallback ret) throws Exception {
		BufferedReader br = new BufferedReader(new FileReader(f));
		while (true) {
			String l = br.readLine();
			if (l == null)
				break;
			if (l.trim().isEmpty())
				continue;
			String[] parts = tabDiv.split(l);
			ret.proteinPair(parts[0], parts[1], Double.parseDouble(parts[2]), Integer.parseInt(parts[3]), 
					Integer.parseInt(parts[4]), Integer.parseInt(parts[5]), Integer.parseInt(parts[6]), 
					Integer.parseInt(parts[7]), Integer.parseInt(parts[8]), Integer.parseInt(parts[9]), 
					parts[10], Double.parseDouble(parts[11]));
		}
		br.close();
	}

	public static interface ResultCallback {
		public void proteinPair(String qname, String tname, double ident, int alnLen, int mismatch,
				int gapopens, int qstart, int qend, int tstart, int tend, String eval, double bitScore);
	}
}