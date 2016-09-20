package genomeproteomecomparison;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Map;

import us.kbase.common.service.Tuple11;

public class Utils {
    public static int compare(long x, long y) {
        return (x < y) ? -1 : ((x == y) ? 0 : 1);
    }

    public static int compare(int x, int y) {
        return (x < y) ? -1 : ((x == y) ? 0 : 1);
    }
    
    public static long copy(InputStream from, OutputStream to) throws IOException {
    	byte[] buf = new byte[10000];
    	long total = 0;
    	while (true) {
    		int r = from.read(buf);
    		if (r == -1) {
    			break;
    		}
    		to.write(buf, 0, r);
    		total += r;
    	}
    	return total;
    }
    
    public static String getRefFromObjectInfo(Tuple11<Long, String, String, String, 
            Long, String, Long, String, String, Long, Map<String,String>> info) {
        return info.getE7() + "/" + info.getE1() + "/" + info.getE5();
    }
}
