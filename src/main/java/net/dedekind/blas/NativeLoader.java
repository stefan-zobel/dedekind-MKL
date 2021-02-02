/*
 * Copyright 2021 Stefan Zobel
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package net.dedekind.blas;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;
import java.nio.channels.ReadableByteChannel;
import java.util.LinkedHashSet;

/**
 * This class is used to load the dedekind-mkl shared library either from the
 * file system or from within the jar.
 */
final class NativeLoader {

    private static final String NATIVE_NAME = "dedekind-mkl";
    private static final String OS = System.getProperty("os.name", "").toLowerCase();
    private static final String LIB_PATH = System.getProperty("java.library.path", "");
    private static final String TMP_DIR = System.getProperty("java.io.tmpdir", "");
    private static final String PROP_DIR = System.getProperty("dedekind.sharedlib.dir", "");
    private static final String ENV_DIR = System.getenv("DEDEKIND_SHAREDLIB_DIR");
    private static final LinkedHashSet<String> ALL_DIRS;
    private static boolean loaded = false;
    private static String location = null;

    static {
        LinkedHashSet<String> directories = new LinkedHashSet<>();
        if (PROP_DIR != null && !PROP_DIR.isEmpty()) {
            String dir = checkDirectory(PROP_DIR);
            if (dir != null) {
                directories.add(dir);
            }
        }
        if (ENV_DIR != null && !ENV_DIR.isEmpty()) {
            String dir = checkDirectory(ENV_DIR);
            if (dir != null) {
                directories.add(dir);
            }
        }
        if (TMP_DIR != null && !TMP_DIR.isEmpty()) {
            String dir = checkDirectory(TMP_DIR);
            if (dir != null) {
                directories.add(dir);
            }
        }
        if (LIB_PATH != null && !LIB_PATH.isEmpty()) {
            String javaLibPath[] = LIB_PATH.split(File.pathSeparator);
            if (javaLibPath != null && javaLibPath.length > 0) {
                for (int i = 0; i < javaLibPath.length; ++i) {
                    String dir = checkDirectory(javaLibPath[i]);
                    if (dir != null) {
                        directories.add(dir);
                    }
                }
            }
        }
        ALL_DIRS = directories;
    }

    static synchronized boolean isLoaded() {
        return loaded;
    }

    static synchronized String getLocation() {
        return location;
    }

    static synchronized void load() {
        if (loaded) {
            return;
        }
        String libName = getLibraryFileName();
        // attempt to load from file system
        for (String dir : ALL_DIRS) {
            File file = new File(dir, libName).getAbsoluteFile();
            if (file.exists() && file.isFile() && load(file)) {
                set(file);
                return;
            }
        }
        // extract from jar file
        try {
            File file = extract(libName);
            if (file != null && file.exists() && file.isFile() && load(file)) {
                set(file);
                return;
            }
        } catch (IOException e) {
            throw new ExceptionInInitializerError(e);
        }

        throw new ExceptionInInitializerError("Unable to load " + libName + " from " + ALL_DIRS.toString());
    }

    private static boolean load(File file) {
        try {
            System.load(file.getAbsolutePath());
            return true;
        } catch (UnsatisfiedLinkError e) {
            if (TMP_DIR != null && file.getAbsolutePath().startsWith(TMP_DIR)) {
                // may be corrupted
                try {
                    file.delete();
                } catch (Exception ignore) {
                }
            }
            return false;
        } catch (SecurityException e) {
            return false;
        } catch (Throwable e) {
            throw new ExceptionInInitializerError(e);
        }
    }

    @SuppressWarnings("resource")
    private static File extract(String fileName) throws IOException {
        String filePath = "/" + fileName;
        URL url = NativeLoader.class.getResource(filePath);
        if (url == null) {
            return null;
        }
        File file = createFile(fileName);
        InputStream in = null;
        FileChannel dest = null;
        FileOutputStream fos = null;
        try {
            in = NativeLoader.class.getResourceAsStream(filePath);
            fos = new FileOutputStream(file);
            dest = fos.getChannel();
            ReadableByteChannel src = Channels.newChannel(in);
            dest.transferFrom(src, 0L, 0x7fffffffffffffffL);
            return file;
        } catch (IOException | SecurityException e) {
            return null;
        } finally {
            try {
                if (dest != null) {
                    dest.close();
                }
                if (fos != null) {
                    fos.close();
                }
                if (in != null) {
                    in.close();
                }
            } catch (Exception ignore) {
            }
        }
    }

    private static File createFile(String libName) throws IOException {
        String dir = null;
        if (!ALL_DIRS.isEmpty()) {
            // simply use the first one
            dir = ALL_DIRS.iterator().next();
        }
        if (dir != null) {
            File file = new File(dir, libName);
            if (file.exists() && !file.isFile()) {
                throw new IllegalArgumentException(file.getAbsolutePath() + " is not a file");
            }
            if (!file.exists()) {
                file.createNewFile();
            }
            return file;
        }
        return null;
    }

    private static void set(File file) {
        loaded = true;
        location = file.getAbsolutePath();
    }

    private static String checkDirectory(String dirPath) {
        if (dirPath != null && !dirPath.isEmpty()) {
            File dir = new File(dirPath).getAbsoluteFile();
            if (dir.exists() && dir.isDirectory()) {
                return dir.getAbsolutePath();
            }
        }
        return null;
    }

    private static String getLibraryFileName() {
        return appendLibOsExtension(prependLibOsPrefix(NATIVE_NAME));
    }

    private static String prependLibOsPrefix(String libFileName) {
        if (isLinux()) {
            return "lib" + libFileName;
        }
        return libFileName;
    }

    private static String appendLibOsExtension(String libFileName) {
        if (isLinux()) {
            return libFileName + ".so";
        }
        if (isWindows()) {
            return libFileName + ".dll";
        }
        throw new UnsupportedOperationException();
    }

    private static boolean isWindows() {
        return OS.contains("windows");
    }

    private static boolean isLinux() {
        return OS.contains("linux");
    }
}
