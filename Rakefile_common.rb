#
#
#
%w(colorize fileutils pathname rubygems/package net/http zip zlib uri openssl).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require_relative "./Rakefile_conf.rb"

if (/cygwin|mswin|mingw|bccwin|wince|emx/ =~ RUBY_PLATFORM) != nil then
  #linux
  task :default => [:install_linux]
elsif (/darwin/ =~ RUBY_PLATFORM) != nil then
  #osx
  task :default => [:install_osx]
else
  #windows
  task :default => [:install_windows]
end


def win_vs( bits, year )

  tmp = " -DBITS=#{bits} -DYEAR=#{year} "

  win32_64  = ''
  win32_64_ = '-A Win32'
  case bits
  when /x64/
    win32_64  = ' Win64'
    win32_64_ = ' -A x64'
  end

  case year
  when "2010"
    tmp = 'cmake -G "Visual Studio 10 2010' + win32_64 +'" ' + tmp
  when "2012"
    tmp = 'cmake -G "Visual Studio 11 2012' + win32_64 +'" ' + tmp
  when "2013"
    tmp = 'cmake -G "Visual Studio 12 2013' + win32_64 +'" ' + tmp
  when "2015"
    tmp = 'cmake -G "Visual Studio 14 2015' + win32_64 +'" ' + tmp
  when "2017"
    tmp = 'cmake -G "Visual Studio 15 2017' + win32_64 +'" ' + tmp
  when "2019"
    tmp = 'cmake -G "Visual Studio 16 2019"' + win32_64_ + tmp
  else
    puts "Visual Studio year #{year} not supported!\n";
    return ""
  end
  return tmp
end