#!/usr/bin/ruby -w

require 'fileutils'

FileUtils.rm_rf   "src"
FileUtils.mkdir_p "src"
FileUtils.mkdir_p "src/Utils/fmt"
FileUtils.mkdir_p "src/Utils/zstream"
FileUtils.mkdir_p "src/Utils/Eigen"
FileUtils.mkdir_p "src/Utils/mingw-std-threads"
FileUtils.mkdir_p "src/tests"

Dir.glob('bin/*.mex*').each { |file| File.delete(file) }

lst = Dir["../src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../src/tests/*.cxx"]
lst.each do |filename|
  FileUtils.cp filename, "./src/tests/" + File.basename(filename);
end

lst = Dir["../submodules/Utils/src/*.cc"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/fmt/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/fmt/*.c*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/fmt/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/zstream/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/zstream/" + File.basename(filename);
end
lst = Dir["../submodules/Utils/src/Utils/mingw-std-threads/*.h*"]
lst.each do |filename|
  FileUtils.cp filename, "./src/Utils/mingw-std-threads/" + File.basename(filename);
end

FileUtils.cp_r "../submodules/Utils/src/Eigen", "src"
FileUtils.cp "../license.txt", "license.txt"
FileUtils.cp "../license_3rd.txt", "license_3rd.txt"
